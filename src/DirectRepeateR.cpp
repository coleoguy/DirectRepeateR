#include <Rcpp.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <limits>
#include <string_view>
#include <cctype>
#include <cstring>
#include <cstdint>

#ifdef __GLIBC__
#include <malloc.h>  // for malloc_trim() if using GNU C library
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// A small struct for storing each repeat row within a chromosome.
// -----------------------------------------------------------------------------
struct RepeatRow {
  int startPos;
  int endPos;
  int matchPos;
  int matchEndPos;
  int nextSp;
  int nextMp;
};

// -----------------------------------------------------------------------------
// A key (sp, mp) for adjacency lookups in a hash map.
// -----------------------------------------------------------------------------
struct Key {
  int sp;
  int mp;
  
  bool operator==(const Key &other) const {
    return (sp == other.sp && mp == other.mp);
  }
};

struct KeyHash {
  std::size_t operator()(const Key &k) const {
    // Combine sp, mp. Shift + XOR to reduce collisions.
    return (std::hash<int>()(k.sp) << 1) ^ std::hash<int>()(k.mp);
  }
};

// -----------------------------------------------------------------------------
// Pattern matcher: returns all 1-based match positions of 'pattern'
// in 'text'. For patterns of at least 8 characters, an unaligned 64-bit
// load compares the first 8 bytes in a single instruction (DNA has a
// 4-letter alphabet, so a one-byte prefilter hits every ~4 positions;
// an 8-byte prefilter virtually never false-positives). The remainder
// is verified with memcmp. Preserves the original semantics: 1-based
// offsets, and the scan advances by the pattern length after each match.
// -----------------------------------------------------------------------------
std::vector<int> findMatchesFixed(std::string_view text, std::string_view pattern) {
  std::vector<int> positions;
  if (pattern.empty() || text.empty() || pattern.size() > text.size()) {
    return positions;
  }
  
  const size_t plen = pattern.size();
  const char  *base = text.data();
  const size_t last = text.size() - plen;  // last valid start offset
  
  if (plen >= 8) {
    std::uint64_t pat8;
    std::memcpy(&pat8, pattern.data(), 8);
    const char  *ptail = pattern.data() + 8;
    const size_t tlen  = plen - 8;
    
    for (size_t i = 0; i <= last; ) {
      std::uint64_t txt8;
      std::memcpy(&txt8, base + i, 8);  // safe: i + plen <= text.size()
      if (txt8 == pat8 && std::memcmp(base + i + 8, ptail, tlen) == 0) {
        positions.push_back(static_cast<int>(i + 1)); // 1-based
        i += plen;  // Skip ahead by pattern length on a match
      } else {
        ++i;
      }
    }
  } else {
    for (size_t i = 0; i <= last; ) {
      if (std::memcmp(base + i, pattern.data(), plen) == 0) {
        positions.push_back(static_cast<int>(i + 1)); // 1-based
        i += plen;
      } else {
        ++i;
      }
    }
  }
  return positions;
}

// -----------------------------------------------------------------------------
// Union-Find (Disjoint Set) helpers
// -----------------------------------------------------------------------------
// Iterative find with path compression (recursion could overflow the
// stack on multi-Mb tandem arrays).
int findRoot(std::vector<int> &parent, int x) {
  int root = x;
  while (parent[root] != root) {
    root = parent[root];
  }
  while (parent[x] != root) {
    int next = parent[x];
    parent[x] = root;
    x = next;
  }
  return root;
}

void unionSets(std::vector<int> &parent, std::vector<int> &rank, int x, int y) {
  int rx = findRoot(parent, x);
  int ry = findRoot(parent, y);
  if (rx == ry) return;
  
  if (rank[rx] > rank[ry]) {
    parent[ry] = rx;
  } else if (rank[rx] < rank[ry]) {
    parent[rx] = ry;
  } else {
    parent[ry] = rx;
    rank[rx]++;
  }
}

// -----------------------------------------------------------------------------
// Process a single chromosome
// -----------------------------------------------------------------------------
void processSingleChrom(const std::string &chromName,
                        const std::string &sequence,
                        int query_length,
                        int maxdist,
                        const std::string &outdir)
{
  // Guard against int overflow on chromosomes > 2^31 - 1 bases
  if (sequence.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
    stop("Chromosome " + chromName +
         " is longer than 2^31 - 1 bases, which is not currently supported.");
  }
  int seq_length = static_cast<int>(sequence.size());
  if (seq_length < query_length) {
    Rcout << "Chromosome " << chromName 
          << " is shorter than query_length; skipping\n";
    return;
  }
  
  // Vector to store found repeats
  std::vector<RepeatRow> repeats;
  repeats.reserve(seq_length / query_length + 100);
  
  int num_chunks = seq_length / query_length;
  
  for (int c = 0; c < num_chunks; ++c) {
    int start_pos = c * query_length + 1; 
    int end_pos   = start_pos + query_length - 1;
    if (end_pos > seq_length) break;
    
    int upper_bound = std::min(end_pos + maxdist, seq_length);
    int search_len  = upper_bound - end_pos;
    if (search_len <= 0) continue;
    
    std::string_view pattern(sequence.data() + (start_pos - 1), query_length);
    
    // Skip chunks containing N (assembly gaps): N-runs would otherwise
    // match each other and generate large false repeat blocks.
    if (pattern.find('N') != std::string_view::npos) continue;
    
    std::string_view search_field(sequence.data() + end_pos, search_len);
    
    std::vector<int> matches = findMatchesFixed(search_field, pattern);
    for (int offset : matches) {
      int abs_match_pos = end_pos + offset; 
      int match_end_pos = abs_match_pos + query_length - 1;
      
      RepeatRow rr;
      rr.startPos    = start_pos;
      rr.endPos      = end_pos;
      rr.matchPos    = abs_match_pos;
      rr.matchEndPos = match_end_pos;
      rr.nextSp      = end_pos + 1;
      rr.nextMp      = match_end_pos + 1;
      repeats.push_back(rr);
    }
  }
  
  if (repeats.empty()) {
    Rcout << "No repeats found for chromosome " << chromName << "\n";
    return;
  }
  
  // Build hash map for adjacency
  std::unordered_map<Key, int, KeyHash> rowIndex;
  rowIndex.reserve(repeats.size());
  
  for (int i = 0; i < static_cast<int>(repeats.size()); ++i) {
    Key k{repeats[i].startPos, repeats[i].matchPos};
    rowIndex[k] = i;
  }
  
  // Union-find prep
  int n = repeats.size();
  std::vector<int> parent(n), rank(n, 0);
  for (int i = 0; i < n; i++) {
    parent[i] = i;
  }
  
  // Adjacency union
  for (int i = 0; i < n; ++i) {
    Key adj{repeats[i].nextSp, repeats[i].nextMp};
    auto it = rowIndex.find(adj);
    if (it != rowIndex.end()) {
      unionSets(parent, rank, i, it->second);
    }
  }
  
  // Path compression
  for (int i = 0; i < n; i++) {
    parent[i] = findRoot(parent, i);
  }
  
  // Condense connected components
  struct Condensed {
    int startPos;
    int endPos;
    int matchPos;
    int matchEndPos;
  };
  
  std::unordered_map<int, Condensed> condensed_map;
  condensed_map.reserve(n / 2);
  
  for (int i = 0; i < n; ++i) {
    int comp = parent[i];
    if (condensed_map.find(comp) == condensed_map.end()) {
      condensed_map[comp] = Condensed {
        repeats[i].startPos,
        repeats[i].endPos,
        repeats[i].matchPos,
        repeats[i].matchEndPos
      };
    } else {
      auto &c = condensed_map[comp];
      c.startPos    = std::min(c.startPos, repeats[i].startPos);
      c.endPos      = std::max(c.endPos, repeats[i].endPos);
      c.matchPos    = std::min(c.matchPos, repeats[i].matchPos);
      c.matchEndPos = std::max(c.matchEndPos, repeats[i].matchEndPos);
    }
  }
  
  // Write CSV for this chromosome
  std::string out_file = outdir + "/" + chromName + "_condensed.csv";
  std::ofstream outfile(out_file);
  if (!outfile.is_open()) {
    stop("Cannot open output file: " + out_file);
  }
  
  // Write header
  outfile << "Start,End,Match_Start,Match_End\n";
  for (const auto &kv : condensed_map) {
    const Condensed &row = kv.second;
    outfile << row.startPos    << ","
            << row.endPos      << ","
            << row.matchPos    << ","
            << row.matchEndPos << "\n";
  }
  outfile.close();
  
  Rcout << "Chromosome " << chromName 
        << " condensed results -> " << out_file << "\n";
  
  // Release memory
  std::vector<RepeatRow>().swap(repeats);
  std::unordered_map<Key, int, KeyHash>().swap(rowIndex);
  std::vector<int>().swap(parent);
  std::vector<int>().swap(rank);
  std::unordered_map<int, Condensed>().swap(condensed_map);
  
#ifdef __GLIBC__
  malloc_trim(0);
#endif
}

// -----------------------------------------------------------------------------
// Sanitize a chromosome name so it can be used as a file name.
// Replaces any character other than [A-Za-z0-9._-] with '_'.
// (NCBI-style headers like ">gi|...|ref|NC_003279.8|" would otherwise
// produce invalid file paths and be silently dropped.)
// -----------------------------------------------------------------------------
std::string sanitizeChromName(const std::string &name) {
  std::string out = name;
  for (char &c : out) {
    bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
              (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
    if (!ok) c = '_';
  }
  return out;
}

// -----------------------------------------------------------------------------
// Process the FASTA by chromosome
// -----------------------------------------------------------------------------
void processFastaByChrom(const std::string &fasta_path,
                         int query_length,
                         int maxdist,
                         const std::string &outdir)
{
  std::ifstream infile(fasta_path);
  if (!infile.is_open()) {
    stop("Cannot open FASTA file: " + fasta_path);
  }
  
  std::string line;
  std::string current_chrom;
  std::string seqbuf;
  
  while (std::getline(infile, line)) {
    // Strip trailing carriage return from Windows (CRLF) files: it would
    // otherwise corrupt matching and shift coordinates.
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) continue;
    
    if (line[0] == '>') {
      if (!current_chrom.empty()) {
        processSingleChrom(current_chrom, seqbuf, query_length, maxdist, outdir);
        seqbuf.clear();
        seqbuf.shrink_to_fit();
      }
      // Extract name from header
      std::string header = line.substr(1);
      size_t spacePos = header.find_first_of(" \t");
      if (spacePos != std::string::npos) {
        header = header.substr(0, spacePos);
      }
      current_chrom = sanitizeChromName(header);
    } else {
      // Uppercase on read-in so soft-masked (lowercase) genomes are
      // matched case-insensitively.
      std::transform(line.begin(), line.end(), line.begin(),
                     [](unsigned char ch) { return std::toupper(ch); });
      seqbuf += line;
    }
  }
  
  // final chromosome
  if (!current_chrom.empty()) {
    processSingleChrom(current_chrom, seqbuf, query_length, maxdist, outdir);
    seqbuf.clear();
    seqbuf.shrink_to_fit();
  }
  
  infile.close();
  
#ifdef __GLIBC__
  malloc_trim(0);
#endif
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
void run_combined_cpp(const std::string &fasta_path,
                      int query_length,
                      int maxdist,
                      const std::string &outdir)
{
  // Create the intermediate-results folder if it does not exist.
  Function dirCreate("dir.create");
  dirCreate(outdir, _["showWarnings"] = false, _["recursive"] = true);
  
  processFastaByChrom(fasta_path, query_length, maxdist, outdir);
  Rcout << "All chromosomes processed. CSVs placed in '" << outdir << "'\n";
}
