#ifndef __CHROM_LOCI_H
#define __CHROM_LOCI_H

#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <algorithm>

#include "Error.h"
#include "hts_utils.h"

// A single genomic int32_terval
class posLocus {
 public:
  int32_t beg1; // includes 1-based, excludes 0-based
  int32_t end0; // excludes 1-based, includes 0-based

  // constructor
  // b : beg1 value
  // e : end0 value
  posLocus(int32_t b, int32_t e) : beg1(b), end0(e) {}

  // compare between genomeLocus
  // l : rhs argument of the same type
  // returns true iff this < l
  bool operator< (const posLocus& l) const {
    if ( beg1 == l.beg1 ) { return end0 < l.end0; }
    else return beg1 < l.beg1;    
  }

  // returns length of the interval
  unsigned long length() const { return end0-beg1+1; }

  // returns the total overlapping intervals
  int32_t overlapBases(int32_t _beg1, int32_t _end0) const {
    if ( ( beg1 <= end0 ) && ( _beg1 <= end0 ) ) {
      return ( _beg1 < end0 ? _beg1 : end0 ) - ( beg1 < end0 ? end0 : beg1 ) + 1;
    }
    else return 0;    
  }

  int32_t overlapBases(const posLocus& l) const {
    return overlapBases(l.beg1, l.end0);
  }

  bool overlaps(int32_t _beg1, int32_t _end0) const {
    if ( ( beg1 <= _end0 )  && ( _beg1 <= end0 ) )
      return true;
    else 
      return false;
  }

  // check overlap with other locus
  bool overlaps (const posLocus& l) const {
    return overlaps(l.beg1, l.end0);
  }

  // merge two locus if possible
  bool merge (const posLocus& l) {
    if ( ( beg1-1 <= l.end0 )  && ( l.beg1-1 <= end0 ) ) {
      if ( l.beg1 < beg1 ) beg1 = l.beg1;
      if ( l.end0 > end0 ) end0 = l.end0;
      return true;
    }
    else {
      return false;
    }
  }

  // check whether the interval contains a particular position in 1-based coordinate
  bool contains1(int32_t pos1 = INT_MAX) const {
    return ( ( pos1 >= beg1 ) && ( pos1 <= end0 ) );
  }

  // check whether the interval contains a particular position in 0-based coordinate  
  bool contains0(int32_t pos0) const { return contains1(pos0+1); }

  // parse a string in [chr]:[beg1]-[end0] format 
  static bool parseRegion(const char* region, std::string& chrom, int32_t& beg1, int32_t& end0) {
    char buf[255];
    strcpy(buf,region);
    const char* pcolon = strchr(region,':');
    const char* pminus = strchr(pcolon+1,'-');
    if ( pcolon == NULL )
      error("Cannot parse %s in genomeLocus::genomeLocus()");
    chrom.assign(region,0,pcolon-region);
    beg1 = atoi(pcolon+1);
    if ( pminus == NULL ) end0 = INT_MAX;
    else {
      end0 = atoi(pminus+1);
      if ( end0 == 0 ) end0 = INT_MAX;
    }
  }
};

// posLoci represents a set of posLocus objects residing in the same chromosome
class posLoci {
public:
  std::set<posLocus> loci;
  bool overlapResolved;

  posLoci() : overlapResolved(true) {}
  
  // reset the iterator to the beginning position
  inline void rewind(std::set<posLocus>::iterator& it) {
    it = loci.begin();
  }

  // check whether the iterator reached at the end
  inline bool isend(std::set<posLocus>::iterator& it) {
    return it == loci.end();
  }

  // advance the iterator to the next element
  inline bool next(std::set<posLocus>::iterator& it) {
    if ( it == loci.end() ) return false;
    return ( it++ != loci.end() );
  }

  // check whether the loci is empty
  inline bool empty() { return loci.empty(); }

  // add a new posLocus object
  bool add(int32_t beg1, int32_t end0) {
    std::pair<std::set<posLocus>::iterator, bool> ret = loci.insert(posLocus(beg1,end0));
    overlapResolved = false;
    return ret.second;
  }

  // resolves overlapping intervals,
  // and returns the total number of unique intervals
  int32_t resolveOverlaps(int32_t& maxLength) {
    std::set<posLocus>::iterator prev, it;
    for(it = loci.begin(); it != loci.end(); ++it) {
      if ( it == loci.begin() ) {
	prev = it;
      }
      else {
	if ( prev->overlaps(*it) ) {
	  posLocus locus = *prev;
	  locus.merge(*it);
	  if ( (int32_t)locus.length() > maxLength ) maxLength = locus.length();
	  loci.erase(it);
	  loci.erase(prev);
	  prev = it = loci.insert(locus).first;
	  ++numMerged;
	}
	else {
	  prev = it;
	}
      }
    }
    overlapResolved = true;
    return numMerged;
  }

  // returns the total length of the intervals
  // if overlaps are not resolved, it will double count redundant regions
  uint64_t totalLength() const {
    uint64_t sz = 0;
    for(it = loci.begin(); it != loci.end(); ++it) {
      sz += it->length();
    }
    return sz;
  }

  // move the current iterator to point the first interval that contains pos1
  bool moveTo(int32_t pos1, std::set<posLocus>::iterator& it) {
    it = loci.lower_bound(posLocus(pos1,pos1));
    return ( it == loci.end() );
  }

  bool contains1(int32_t pos1) {
    std::set<posLocus>::iterator it = loci.lower_bound(posLocus(pos1,pos1));
    it = 
  }
};

// chromLoci is alternative to genomeLoci
class chromLoci {
public:
  std::map<std::string, posLoci> chrom2loci;
  std::map<std::string, posLoci>::iterator chromIt;
  std::set<posLocus>::iterator posIt;
  bool overlapResolved;
  int32_t maxLength;

  chromLoci() : overlapResolved(false), maxLength(0) {}
  
  chromLoci(const char* reg) : overlapResolved(false), maxLength(0) {
    add(reg);
    //resolveOverlaps();
  }

  inline bool isend() { return chromIt == chrom2loci.end(); }  
  inline void rewind() {
    chromIt = chrom2loci.begin();
    while ( !isend() ) {
      if ( chromIt->second.empty() ) { // skip empty chromosome
	++chromIt;
      }
    }
    if ( !isend() )  // contains at least one element
      chromIt->second.rewind(posIt); // points first non-empty elem
  }

  bool next() {
    // check if current chrom is finished iteration
    bool chromChanged = false;
    while( !isend() ) {
      if ( chromIt->second.empty() ) { // skip empty chromosome
	++chromIt;
	chromChanged = true;
      }
      else { // non-empty chromosome
	//if ( chromIt->second.isend(posIt) ) { // this should never happen
	//  error("E[%s:%d %s]: Observed impossible condition", __FILE__, __LINE__, __PRETTY_FUNCTION__);
	//}
	if ( chromChanged ) 
	  postIt->second.rewind(posIt);
	else if ( posIt->second.next(it) ) return true;
	else { // current chromosome is fully iterated
	  ++chromIt;
	  chromChanged = true;	  
	}
      }
    }
    return false; // reached at the end of the chromosome
  }

  bool clear() {
    if ( chrom2loci.empty() ) return false;
    chrom2loci.clear();
    rewind();
    overlapResolved = false;
    maxLength = 0;
    return true;
  }

  int32_t numLocus() const {
    int32_t ret = 0;
    for(std::map<std::string, posLoci>::iterator it = chrom2loci.begin();
	it != chrom2loci.end(); ++it) {
      ret += (int32_t) it->size();
    }
    return ret;
  }

  uint64_t totalLength() const {
    uint64_t sz = 0;
    for(std::map<std::string, posLoci>::iterator it = chrom2loci.begin();
	it != chrom2loci.end(); ++it) {
      sz += it->totalLength();
    }
    return sz;    
  }

  bool hasChrom(const char* chr) {
    return ( chrom2loci.find(chr) != chrom2loci.end() );
  }

  bool moveTo(const char* chr = NULL, int32_t pos1 = INT_MAX ) {
    //notice("[%s:%d %s] (%s, %d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr == NULL ? "NULL" : chr, pos1);

    // do nothing if empty
    if ( chrom2loci.empty() ) return false;

    // see if chromosome needs to be updated
    bool chromUpdated = false;
    if ( ( chromIt == chrom2loci.end() ) || ( chromIt->first.compare(chr) == 0 ) ) { // chromosome was changed
      chromIt = chrom2loci.lower_bound(chr);
      if ( chromIt == chrom2loci.end() ) { // no more elements found
	return false;
      }
      chromUpdated = true;
      if ( chromIt->first.compare(chr) == 0 ) { // chromosome found
	if ( chromIt->second.moveTo(pos1, posIt) ) return true;
	else return next();
      }
    }
  }

  bool add(const char* chr, int32_t beg1, int32_t end0) {
    overlapResolved = false;
    if ( end0-beg1+1 > maxLength ) maxLength = end0-beg1+1;

    return chrom2loci[chr].add(beg1, end0);
  }

  bool add(const char* region) {
    std::string chrom;
    int32_t beg1, end0;
    posLocus::parseRegion(region, chrom, beg1, end0);
    return add(chrom.c_str(), beg1, end0);
  }

  // Resolve overlapping int32_tervals
  int32_t resolveOverlaps() {
    if ( !overlapResolved ) {
      int32_t numMerged = 0;
      for(std::map<std::string, posLoci>::iterator it = chrom2loci.begin();
	  it != chrom2loci.end(); ++it) {
	numMerged += it->resolveOverlaps(maxLnegth);
      }
      return numMerged;
    }
    else return 0;
  }

  bool contains1(const char* chr, int32_t pos1) {
    if ( chrom2loci.empty() ) return false;
    std::map<std::string, posLoci>::iterator it = chrom2loci.find(chr);

    if ( it == chrom2loci.end() ) return false;
    else return it->second.contains1(pos1);
  }  

  bool openBED(const char* file) {
    clear();
    
    htsFile* fp = hts_open(file, "r");
    if ( fp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, file);
    
    kstring_t str = {0,0,0};
    int32_t lstr = 0;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    int32_t i;
    // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
    for( i=0; ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, file);      
      // typically bed files have beg0 and end0
      add(&str.s[fields[0]], atoi(&str.s[fields[1]])+1, atoi(&str.s[fields[2]]));
      //if ( i % verbose == 0 )
      //notice("Processing %d lines to mask at %s:%d-%d in %s..", i, &str.s[fields[0]], &str.s[fields[1]], &str.s[fields[2]], file);
    }
    hts_close(fp);
    
    notice("Processed %d lines from %s, maxLength = %d", i, file, maxLength);
    
    resolveOverlaps();

    notice("After removing overlaps, %d intervals remained, maxLength = %d, total length = %u", (int32_t)loci.size(), maxLength, totalLength());

    rewind();
    return i > 0;
  }
}

// For single chromosome, multiple posLocus elements exists
// Useful for representing exons, for example
class singleChromLoci {
public:
  std::string chrom;
  std::set<posLocus> loci;
};

class multiChromLoci {
public:
  std;:map<std::string, posLoci> chrom2loci;
};

class chromLoci {
 public:
  std::set<std::string> chroms;
  std::set<genomeLocus> loci;
  std::set<genomeLocus>::iterator it;
  bool overlapResolved;
  int32_t maxLength;

  genomeLoci() : overlapResolved(false), maxLength(0) {}
  genomeLoci(const char* reg) : overlapResolved(false), maxLength(0) {
    add(reg);
    resolveOverlaps();
  }

  // functions for iterating each locus
  inline void rewind() { it = loci.begin(); }
  inline bool next() { ++it; return ( it != loci.end() ); }
  inline bool isend() { return ( it == loci.end() );  }
  inline const genomeLocus& currentLocus() { return (*it); }

  // check the size 
  bool empty() { return loci.empty(); }
  
  bool clear() {
    if ( loci.empty() ) return false;
    
    loci.clear();
    it = loci.begin();
    overlapResolved = false;
    maxLength = 0;
    return true;
  }
 
  int32_t numLocus() const { return (int32_t)loci.size(); }

  bool openBED(const char* file) {
    clear();
    
    htsFile* fp = hts_open(file, "r");
    if ( fp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, file);
    
    kstring_t str = {0,0,0};
    int32_t lstr = 0;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    int32_t i;
    // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
    for( i=0; ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, file);      
      // typically bed files have beg0 and end0
      add(&str.s[fields[0]], atoi(&str.s[fields[1]])+1, atoi(&str.s[fields[2]]));
      //if ( i % verbose == 0 )
      //notice("Processing %d lines to mask at %s:%d-%d in %s..", i, &str.s[fields[0]], &str.s[fields[1]], &str.s[fields[2]], file);
    }
    hts_close(fp);
    
    notice("Processed %d lines from %s, maxLength = %d", i, file, maxLength);
    
    resolveOverlaps();

    notice("After removing overlaps, %d intervals remained, maxLength = %d, total length = %u", (int32_t)loci.size(), maxLength, totalLength());

    rewind();

    return i > 0;
  }

  // add a locus
  bool add(const char* chr, int32_t beg1, int32_t end0) {
    overlapResolved = false;
    if ( end0-beg1+1 > maxLength ) maxLength = end0-beg1+1;
    std::pair<std::set<genomeLocus>::iterator, bool> ret = loci.insert(genomeLocus(chr,beg1,end0));
    it = ret.first;
    if ( ret.second )
      chroms.insert(ret.first->chrom);    
    return ret.second;
  }

  // add a locus
  bool add(const char* region) {
    overlapResolved = false;
    std::pair<std::set<genomeLocus>::iterator, bool> ret = loci.insert(genomeLocus(region));
    if ( ret.second )
      chroms.insert(ret.first->chrom);
    it = ret.first;
    int32_t l = ret.first->end0 - ret.first->beg1 + 1;
    if ( l > maxLength ) maxLength = l;
    return ret.second;
  }

  // Resolve overlapping int32_tervals
  int32_t resolveOverlaps() {
    if ( !overlapResolved ) {
      std::set<genomeLocus>::iterator it;
      std::set<genomeLocus>::iterator prev;
      int32_t numMerged = 0;
      for(it = loci.begin(); it != loci.end(); ++it) {
	if ( it != loci.begin() ) {
	  if ( prev->overlaps(*it) ) {
	    // if overlaps, erase both and insert merged one
	    genomeLocus locus = *prev;
	    locus.merge(*it);
	    if ( (int32_t)locus.length() > maxLength ) maxLength = locus.length();
	    loci.erase(it);
	    loci.erase(prev);
	    prev = it = loci.insert(locus).first;
	    ++numMerged;
	  }
	  else {
	    prev = it;
	  }
	}
	else {
	  prev = it;
	}
      }
      overlapResolved = true;
      return numMerged;
    }
    else {
      return 0;
    }
    return 0;
  }

  unsigned long totalLength() const {
    //resolveOverlaps();
    unsigned long sz = 0;
    std::set<genomeLocus>::iterator it2;
    for(it2 = loci.begin(); it2 != loci.end(); ++it2) {
      sz += it2->length();
    }
    return sz;
  }

  bool hasChrom(const char* chr) {
    return ( chroms.find(chr) != chroms.end() );
  }

  bool moveTo(const char* chr = NULL, int32_t pos1 = INT_MAX) {
    notice("[%s:%d %s] (%s, %d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr == NULL ? "NULL" : chr, pos1);
    
    if ( loci.empty() ) return false;

    if ( ( it != loci.end() ) && it->contains1(chr, pos1) ) return true;
    
    if ( chr == NULL ) chr = it->chrom.c_str();
    
    genomeLocus locus(chr, pos1, pos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      notice("beg");
      return (it->contains1(chr,pos1));
    }
    else if ( it == loci.end() ) {
      notice("end");      
      std::set<genomeLocus>::iterator i = it;
      --i;
      if ( i->contains1(chr,pos1) ) { it = i; return true; }
      else { rewind(); return false; }
    }
    else {
      notice("mid");                  
      if ( it->contains1(chr,pos1) ) return true;
      else {
	std::set<genomeLocus>::iterator i = it;
	--i;
	if ( i->contains1(chr,pos1) ) { it = i; return true; }
	else { rewind(); return false; }
      }
    }
  }

  bool contains1(const char* chr, int32_t pos1) {
    if ( loci.empty() ) return false;
    notice("contains1(%s,%d) called", chr, pos1);    
    genomeLocus locus(chr, pos1, pos1);
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= pos1 ) ) {
      if ( it2->end0 >= pos1 ) return true;
      ++it2;
    }
    return false;
  }

  bool overlaps(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;
    
    genomeLocus locus(chr, overlapResolved ? beg1 : beg1-maxLength, overlapResolved ? beg1 : beg1-maxLength);
    if ( loci.empty() ) return false;
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= end0 ) ) {
      if ( ( it2->beg1 <= end0 ) && ( beg1 <= it2->end0 ) )
	return true;
      ++it2;
    }
    //notice("%s:%d-%d",it2->chrom.c_str(),it2->beg1,it2->end0);
    return false;
  }

  bool contains(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;

    resolveOverlaps();
    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= end0 ) ) {
      if ( ( it2->beg1 <= beg1 ) && ( end0 <= it2->end0 ) )
	return true;
      ++it2;
    }
    return false;    
  }
};

// Collection of genomic locus
template <class T>
class genomeLocusMap {
 public:
  std::set<std::string> chroms;  
  std::map<genomeLocus,T> loci;
  typename std::map<genomeLocus,T>::iterator it;
  int32_t maxLength;

 genomeLocusMap() : maxLength(0) { it = loci.end(); }
  genomeLocusMap(const char* reg, const T& val) : maxLength(0) {
    add(reg, val);
    it = loci.end();
  }

  // functions for iterating each locus
  void rewind() { it = loci.begin(); }
  bool next() { ++it; return ( it != loci.end() ); }
  bool isend() { return ( it == loci.end() ); }
  const genomeLocus& currentLocus() { return (*it); }

  // check the size 
  bool empty() { return loci.empty(); }
  
  bool clear() {
    if ( loci.empty() ) return false;
    loci.clear();
    it = loci.begin();
    maxLength = 0;
    return true;
  }
    
  int32_t numLocus() const { return (int32_t)loci.size(); }

  // add a locus
  bool add(const char* chr, int32_t beg1, int32_t end0, const T& val) {
    if ( end0-beg1+1 > maxLength ) maxLength = end0-beg1+1;
    std::pair<typename std::map<genomeLocus,T>::iterator, bool> ret = loci.insert(std::pair<genomeLocus,T>(genomeLocus(chr,beg1,end0),val));
    it = ret.first;
    if ( ret.second )
      chroms.insert(chr);
    return ret.second;
  }
  
  // add a locus
  bool add(const char* region, const T& val) {
    std::pair<typename std::map<genomeLocus,T>::iterator, bool> ret = loci.insert(std::pair<genomeLocus,T>(region,val));
    int32_t l = ret.first->end0 - ret.first->beg1 + 1;
    if ( ret.second )
      chroms.insert(ret.first->chrom);    
    if ( l > maxLength ) maxLength = l;
    return ret.second;
  }

  unsigned long totalLength() const {
    unsigned long sz = 0;
    typename std::map<genomeLocus,T>::iterator it2;
    for(it2 = loci.begin(); it2 != loci.end(); ++it2) {
      sz += it2->first.length();
    }
    return sz;
  }

  bool moveTo(const char* chr = NULL, int32_t pos1 = INT_MAX) {
    notice("[%s:%d %s] (%s, %d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr == NULL ? "NULL" : chr, pos1);

    if ( loci.empty() ) return false;    

    if ( ( it != loci.end() ) && it->first.contains1(chr, pos1) ) return true;    
    
    if ( chr == NULL ) chr = it->first.chrom.c_str();

    genomeLocus locus(chr, pos1, pos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      return (it->first.contains1(chr,pos1));
    }
    else if ( it == loci.end() ) {
      typename std::map<genomeLocus,T>::iterator i = it;
      --i;
      if ( i->first.contains1(chr,pos1) ) { it = i; return true; }
      else { rewind(); return false; }
    }
    else {
      if ( it->first.contains1(chr,pos1) ) return true;
      else {
	typename std::map<genomeLocus,T>::iterator i = it;
	--i;
	if ( i->first.contains1(chr,pos1) ) { it = i; return true; }
	else { rewind(); return false; }
      }
    }
  }

  bool contains1(const char* chr, int32_t pos1) {
    notice("contains1(%s,%d) called", chr, pos1);
    genomeLocus locus(chr, pos1, pos1);
    typename std::map<genomeLocus,T>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->first.chrom == chr ) && ( it2->first.beg1 <= pos1 ) ) {
      if ( it2->first.end0 >= pos1 ) return true;
      ++it2;
    }
    return false;
  }

  bool overlaps(const char* chr, int32_t beg1, int32_t end0) {
    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    if ( loci.empty() ) return false;
    typename std::map<genomeLocus,T>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->first.chrom != chr ) ++it2;
    while( it2 != loci.end() && ( it2->first.chrom == chr ) && ( it2->first.beg1 <= end0 ) ) {
      if ( ( it2->first.beg1 <= end0 ) && ( beg1 <= it2->first.end0 ) )
	return true;
      ++it2;
    }
    return false;
  }

  bool contains(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;

    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    typename std::map<genomeLocus,T>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->first.chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->first.chrom == chr ) && ( it2->first.beg1 <= end0 ) ) {
      if ( ( it2->first.beg1 <= beg1 ) && ( end0 <= it2->first.end0 ) )
	return true;
      ++it2;
    }
    return false;    
  }
};

#endif
