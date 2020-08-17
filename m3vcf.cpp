#include <cassert>
#include <cstring>
#include <algorithm>
#include <map>

#include "m3vcf.h"
#include "Error.h"
#include "hts_utils.h"

void m3vcf_block_t::init(const std::vector<int32_t>& ihaps) {
  // clear up the current block
  ireps.clear();
  sites.clear();
  variants.clear();
  
  assert(ihaps.size() % 2 == 0 ); // only allows diploids for now
  nsamples = (int32_t)ihaps.size() / 2;
  ireps = ihaps;
  //std::copy( ihaps.begin(), ihaps.end(), ireps.begin() );
  int32_t maxrep = *std::max_element( ireps.begin(), ireps.end() );
  nreps = maxrep + 1;
}

int32_t m3vcf_block_t::add_variant(const m3vcf_site_t& site, const std::vector<int32_t>& inrefs) {
  sites.push_back(site);
  variants.push_back(inrefs);
}

bool m3vcf_block_t::write(htsFile* wf) {
  //notice("Entering m3vcf_block_t::write %x", wf, chrom.c_str());  
  // write the block header
  if ( variants.empty() ) return false; // do not write anything if there is no variant
  
  hprintf(wf, "%s\t%d-%d\t<BLOCK>\t.\t.\t.\t.\tVARIANTS=%zu;REPS=%d\t.", chrom.c_str(), sites[0].pos, sites.back().pos, variants.size(), nreps);
  //hprintf(wf, "%s", chrom.c_str());
  //notice("foo");    
  for(int32_t i=0; i < nsamples; ++i) {
    hprintf(wf, "\t%d|%d", ireps[i+i], ireps[i+i+1]);
  }
  hprintf(wf, "\n");

  //notice("Finished writing block header");

  // write each variant
  for(int32_t i=0; i < (int32_t)variants.size(); ++i) {
    const m3vcf_site_t& s = sites[i];
    hprintf(wf, "%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s", s.chrom.c_str(), s.pos, s.varid.c_str(), s.ref.c_str(), s.alt.c_str(), s.qual.c_str(), s.filt.c_str(), s.info.c_str());
    const std::vector<int32_t>& v = variants[i];
    if ( v.size() > 0 ) {
      hprintf(wf, "\t%d", v[0]);
      for(int32_t j=1; j < (int32_t)v.size(); ++j)
	hprintf(wf, ",%d",v[j]);
    }
    hprintf(wf,"\n");
  }
  return true;
}

m3vcf_block_t* m3vcf_block_t::create_subset(const std::vector<int32_t>& isamples, bool removeMono) {
  std::vector<int32_t> new_ireps;
  std::map<int32_t,int32_t> irep_old2new;
  for(int32_t i=0; i < (int32_t)isamples.size(); ++i) {
    int32_t i1 = ireps[isamples[i]*2];
    int32_t i2 = ireps[isamples[i]*2+1];
    std::map<int32_t,int32_t>::iterator it = irep_old2new.find(i1);
    if ( it == irep_old2new.end() ) { // never seen before
      int32_t newid = irep_old2new.size();
      irep_old2new[i1] = newid;
      new_ireps.push_back(newid);
    }
    else {
      new_ireps.push_back(it->second);
    }

    it = irep_old2new.find(i2);
    if ( it == irep_old2new.end() ) { // never seen before
      int32_t newid = irep_old2new.size();
      irep_old2new[i2] = newid;
      new_ireps.push_back(newid);
    }
    else {
      new_ireps.push_back(it->second);
    }    
  }

  // create a new block
  m3vcf_block_t* sblock = new m3vcf_block_t(new_ireps);
  sblock->chrom = chrom;

  //notice("%s %d %d",chrom.c_str(), pos_beg, pos_end);
  //notice("%s %d %d",sblock->chrom.c_str(), sblock->pos_beg, sblock->pos_end);  

  // add variants
  for(int32_t i = 0; i < (int32_t)variants.size(); ++i) {
    std::vector<int32_t> new_alleles(new_ireps.size(),0);
    int32_t idx = 0;
    for(int32_t j=0; j < (int32_t)variants[i].size(); ++j) {
      idx += variants[i][j];
      std::map<int32_t,int32_t>::iterator it = irep_old2new.find(idx);
      if ( it != irep_old2new.end() )
	new_alleles[it->second] = 1;
    }
    std::vector<int32_t> new_inrefs;
    int32_t prev = 0;
    for(int32_t j=0; j < (int32_t)new_alleles.size(); ++j) {
      if ( new_alleles[j] > 0 ) {
	new_inrefs.push_back(j - prev);
	prev = j;
      }
    }
    //notice("i = %d, sites.size() = %zu, new_inrefs.size() = %zu", i, sites.size(), new_inrefs.size());

    // remove monomorphic unless if it is at the boundary
    if ( ( !new_inrefs.empty() ) || ( !removeMono ) || ( i == 0 ) || ( i == (int32_t)variants.size()-1 ) )
      sblock->add_variant(sites[i], new_inrefs);
  }

  // if a block does not contain any variant, call as unlucky and throws an error
  if ( sblock->sites.empty() )
    error("No variant is left in block %d-%d", pos_beg, pos_end);

  // determine boundary positions of the new blocks
  sblock->pos_beg = sblock->sites.front().pos;
  sblock->pos_end = sblock->sites.back().pos;

  if ( ( sblock->pos_beg != pos_beg ) || ( sblock->pos_end != pos_end ) ) { // if the subsetted block's boundary position has changed
    if ( !removeMono )
      error("Something went wrong. Missing variants even though --remove-mono = false. Original block was %d-%d and the new block is %d-%d", pos_beg, pos_end, sblock->pos_beg, sblock->pos_end);
    else
      notice("WARNING: Removed variants at boundary of block %d-%d", pos_beg, pos_end);
  }

  return sblock;
}

int32_t m3vcf_reader_t::open(const char* filename) {
  if ( !tr.open(filename) ) {
    error("Cannot open file %s for reading", filename); 
  }

  // read headers
  while( tr.read_line() > 0 ) {
    const char* token1 = tr.str_field_at(0);
    if ( token1[0] == '#' ) {
      if ( token1[1] == '#' ) {
	// this is a meta-line, reconstruct header
	std::string hdr_line(token1);
	for(int32_t i=1; i < tr.nfields; ++i) {
	  hdr_line += " ";
	  hdr_line += tr.str_field_at(i);
	}
	m3vcf.hdrs.push_back(hdr_line);
      }
      else {
	// this is a header line, and certain names are expected, skip checking for now..
	for(int32_t i=9; i < tr.nfields; ++i) {
	  m3vcf.inds.push_back(tr.str_field_at(i)); // add sample IDs;
	}
	break;
      }
    }
    else {
      error("Observed non-header token %s in file %s when a header/meta line is expected", token1, filename);    
    }
  }

  return (int32_t)m3vcf.inds.size();
}

int32_t m3vcf_reader_t::read_block() {
  // create a new block if needed
  if ( m3vcf.blocks.empty() || m3vcf.keep_in_memory ) {
    m3vcf_block_t* new_block = new m3vcf_block_t;
    m3vcf.blocks.push_back(new_block);
  }

  m3vcf_block_t* p_block = m3vcf.blocks.back();

  // process the first line, we expect to see <BLOCK> with VARIANTS=***;REPS=***
  tr.read_line();
  if ( tr.nfields == 0 )  // EOF was reached
    return 0;
  
  //notice("%d", tr.nfields);

  p_block->chrom.assign(tr.str_field_at(0));
  const char* pos_beg_end = tr.str_field_at(1);
  const char* p_delim = strchr(pos_beg_end, '-');
  p_block->pos_beg = atoi(pos_beg_end);
  p_block->pos_end = atoi(p_delim+1);
  const char* p_info = strchr(tr.str_field_at(7), '=');
  int32_t n_vars = atoi(p_info+1); // Read VARIANTS=***;
  p_info = strchr(p_info+1, '=');
  int32_t n_reps = atoi(p_info+1); // Read REPS=***

  //notice("Reading block %s:%d-%d with %d variants and %d reps", p_block->chrom.c_str(), p_block->pos_beg, p_block->pos_end, n_vars, n_reps);

  std::vector<int32_t> ihaps;
  for(int32_t i=9; i < tr.nfields; ++i) {
    ihaps.push_back(atoi(tr.str_field_at(i)));
    p_delim = strchr(tr.str_field_at(i), '|');
    ihaps.push_back(atoi(p_delim+1));    
  }
  p_block->init(ihaps);

  // process the block, by reading n_vars variants
  std::vector<int32_t> inrefs;
  for(int32_t i=0; i < n_vars; ++i) {
    //notice("i = %d",i);
    tr.read_line();
    m3vcf_site_t site(tr.str_field_at(0), tr.int_field_at(1), tr.str_field_at(2), tr.str_field_at(3), tr.str_field_at(4), tr.str_field_at(5), tr.str_field_at(6), tr.str_field_at(7));
    inrefs.clear();
    const char* pch = tr.nfields > 8 ? tr.str_field_at(8) : NULL;
    while( pch != NULL ) {
      inrefs.push_back(atoi(pch));
      pch = strchr(pch, ',');
      if ( pch != NULL )      
	++pch;
    }
    p_block->add_variant(site, inrefs);
  }

  return n_vars;
}

bool m3vcf_t::write_hdrs(htsFile* wf) {
  for(int32_t i=0; i < (int32_t)hdrs.size(); ++i) 
    hprintf(wf, "%s\n", hdrs[i].c_str());
  
  hprintf(wf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for(int32_t i=0; i < (int32_t)inds.size(); ++i)
    hprintf(wf, "\t%s", inds[i].c_str());
  
  hprintf(wf, "\n");

  //notice("Finished writing headers");

  return true;
}

bool m3vcf_t::write_block(htsFile* wf, int32_t idx) {
  m3vcf_block_t* p = blocks[idx];
  if ( p == NULL )
    error("Cannot find block %d", idx);
  return p->write(wf);
}

bool m3vcf_t::add_block(m3vcf_block_t* new_block) {
  if ( keep_in_memory || blocks.empty() ) { // add the new block
    blocks.push_back(new_block);
    return true;
  }
  else {
    delete blocks.back(); // remove the existing block
    blocks[blocks.size()-1] = new_block;
  }
}
