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
// Naive pattern matcher: returns all 1-based match positions of 'pattern'
// in 'text'. Uses std::string_view to avoid copying large substrings.
// -----------------------------------------------------------------------------
std::vector<int> findMatchesFixed(std::string_view text, std::string_view pattern) {
  std::vector<int> positions;
  if (pattern.empty() || text.empty()) return positions;
  
  size_t plen = pattern.size();
  size_t tlen = text.size();
  
  for (size_t i = 0; i + plen <= tlen; ) {
    if (std::equal(pattern.begin(), pattern.end(), text.begin() + i)) {
      positions.push_back(static_cast<int>(i + 1)); // 1-based
      i += plen;  // Skip ahead by pattern length on a match
    } else {
      ++i;
    }
  }
  return positions;
}

// -----------------------------------------------------------------------------
// Union-Find (Disjoint Set) helpers
// -----------------------------------------------------------------------------
int findRoot(std::vector<int> &parent, int x) {
  if (parent[x] != x) {
    parent[x] = findRoot(parent, parent[x]);
  }
  return parent[x];
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
    Rcout << "Warning: cannot open output file: " << out_file << "\n";
    return;
  }
  
  // Write header
  outfile << "Start_Position,End_Position,Match_Position,Match_End_Position\n";
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
  std::ostringstream seqbuf;
  
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    
    if (line[0] == '>') {
      if (!current_chrom.empty()) {
        std::string chrom_sequence = seqbuf.str();
        processSingleChrom(current_chrom, chrom_sequence, query_length, maxdist, outdir);
        
        // reset
        std::ostringstream empty;
        seqbuf.swap(empty);
      }
      // Extract name from header
      std::string header = line.substr(1);
      size_t spacePos = header.find_first_of(" \t");
      if (spacePos != std::string::npos) {
        header = header.substr(0, spacePos);
      }
      current_chrom = header;
    } else {
      seqbuf << line;
    }
  }
  
  // final chromosome
  if (!current_chrom.empty()) {
    std::string chrom_sequence = seqbuf.str();
    processSingleChrom(current_chrom, chrom_sequence, query_length, maxdist, outdir);
    
    std::ostringstream empty;
    seqbuf.swap(empty);
  }
  
  infile.close();
  
#ifdef __GLIBC__
  malloc_trim(0);
#endif
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
// Hard-coded outdir = "chromosome_results", no user option.
void run_combined_cpp(const std::string &fasta_path,
                      int query_length,
                      int maxdist)
{
  // We always store intermediate CSVs in "chromosome_results".
  std::string outdir = "chromosome_results";
  
  // Create the folder if it does not exist (from C++)
  // We can call R's directory creation function, or just rely on R to do it.
  // We'll attempt to create it here for safety, using system calls or Rcpp 
  // capabilities. But easiest is just to do in R code. 
  // For a small cross-platform approach, let's call R's mkdir:
  Function dirCreate("dir.create");
  dirCreate(outdir, _["showWarnings"] = false, _["recursive"] = true);
  
  processFastaByChrom(fasta_path, query_length, maxdist, outdir);
  Rcout << "All chromosomes processed. CSVs placed in '" << outdir << "'\n";
}
