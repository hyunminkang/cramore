#include "sc_drop_seq.h"
#include "PhredHelper.h"

int32_t sc_dropseq_lib_t::add_snp(int32_t _rid, int32_t _pos, char _ref, char _alt, double _af, double* _gps) {
  snps.resize(nsnps+1);
  snp_umis.resize(nsnps+1);
  sc_snp_t& snp = snps.back();
  snp.rid = _rid;
  snp.pos = _pos;
  snp.ref = _ref;
  snp.alt = _alt;
  snp.af  = _af;
  snp.gps = _gps;

  //if ( _pos % 1000 == 0 ) notice("%lf %lf %lf %lf %lf %lf",_gps[3],_gps[4],_gps[5],_gps[6],_gps[7],_gps[8]);
  
  ++nsnps;
  return nsnps-1;
}

int32_t sc_dropseq_lib_t::add_cell(const char* barcode) {
  std::map<std::string,int32_t>::iterator it = bc_map.find(barcode);
  if ( it == bc_map.end() ) {
    bc_map[barcode] = nbcs;
    bcs.push_back(barcode);
    
    cell_umis.resize( nbcs + 1 );
    cell_totl_reads.resize( nbcs + 1 );
    cell_pass_reads.resize( nbcs + 1 );
    cell_uniq_reads.resize( nbcs + 1 );
    cell_scores.resize(nbcs + 1);
    ++nbcs;
    
    return (nbcs-1);
  }
  else return it->second;
}

bool sc_dropseq_lib_t::add_read(int32_t snpid, int32_t cellid, const char* umi, char allele, char qual) {
  std::map<int32_t,sc_snp_droplet_t*>::iterator it = snp_umis[snpid].find(cellid);
  sc_snp_droplet_t* p_snp_drop = NULL;
  bool ret = false;

  //if ( ( allele > 2 ) || ( qual > 50 ) )
  //  error("umi = %s, allele = %c, qual = %c", umi, allele, qual);

  ++cell_pass_reads[cellid];

  // check if (snp,cell) is empty
  if ( it == snp_umis[snpid].end() ) {
    p_snp_drop = new sc_snp_droplet_t;
    (*p_snp_drop)[umi] = (uint32_t)( ( allele << 24 ) | ( qual << 16 ) | 0x01 );
    snp_umis[snpid][cellid] = p_snp_drop;
    //p_snp_drop = snp_umis[snpid][cellid];
    ret = true;
  }
  // check if (snp,cell,umi) is empty
  else {
    sc_snp_droplet_it_t it2 = it->second->find(umi);
    if ( it2 == it->second->end() ) {
      (*it->second)[umi] = (uint32_t)( ( allele << 24 ) | ( qual << 16 ) | 0x01 );
      ret = true;
    }
    else {
      ++(it->second->at(umi));
    }
    p_snp_drop = it->second;
  }

  std::map<int32_t,sc_snp_droplet_t*>::iterator it3 = cell_umis[cellid].find(snpid);
  if ( it3 == cell_umis[cellid].end() ) {
    cell_umis[cellid][snpid] = p_snp_drop;
  }
  else if ( it3->second != p_snp_drop ) {
    notice("[E:%s] Conflict : Multiple (cellid,snpid) pair %d %d, %x %x", __PRETTY_FUNCTION__, cellid, snpid, it3->second, p_snp_drop);
    notice("%s",p_snp_drop->begin()->first.c_str());
    notice("%d",it3->first);    
    notice("%x",it3->second);
    notice("%u",it3->second->size());
    notice("%s",it3->second->begin()->first.c_str());        
    error("[E:%s] Conflict : Multiple (cellid,snpid) pair %d %d, %s %s", __PRETTY_FUNCTION__, cellid, snpid, p_snp_drop->begin()->first.c_str(), it3->second->begin()->first.c_str());
  }
  if ( ret ) ++cell_uniq_reads[cellid];
  return ret;
}

double calculate_snp_droplet_doublet_GL(sc_snp_droplet_t* ssd, double* gls, double alpha) {
  //double logdenom = 0;
  double tmp;
  gls[0] = gls[1] = gls[2] = gls[3] = gls[4] = gls[5] = gls[6] = gls[7] = gls[8] = 1.0;

  // iterate over each possible bases
  for(sc_snp_droplet_it_t it = ssd->begin(); it != ssd->end(); ++it) {
    uint8_t al = ( it->second >> 24 ) & 0x00ff;
    uint8_t bq = ( it->second >> 16 ) & 0x00ff;

    if ( al == 2 ) continue;

    // 0 : REF / REF -- 1        0
    // 1 : REF / HET -- 1-a/2    a/2
    // 2 : REF / ALT -- 1-a      a
    // 3 : HET / REF -- (1+a)/2  (1-a)/2
    // 4 : HET / HET -- 1/2      1/2
    // 5 : HET / ALT -- (1-a)/2  (1+a)/2
    // 6 : ALT / REF -- a        1-a
    // 7 : ALT / HET -- a/2      1-a/2
    // 8 : ALT / ALT -- 0        1

    gls[0] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1.0           : 0.0          ) + phredConv.phred2Err[bq] / 4. );
    gls[1] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1. - alpha/2. : alpha/2.     ) + phredConv.phred2Err[bq] / 4. );
    gls[2] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1.0 - alpha   : alpha        ) + phredConv.phred2Err[bq] / 4. );
    gls[3] *= ( phredConv.phred2Mat[bq] * (al == 0 ? (1.+alpha)/2. : (1.-alpha)/2.) + phredConv.phred2Err[bq] / 4. );
    gls[4] *= ( phredConv.phred2Mat[bq] * (al == 0 ? .5            : .5           ) + phredConv.phred2Err[bq] / 4. );
    gls[5] *= ( phredConv.phred2Mat[bq] * (al == 0 ? (1.-alpha)/2. : (1.+alpha)/2.) + phredConv.phred2Err[bq] / 4. );
    gls[6] *= ( phredConv.phred2Mat[bq] * (al == 0 ? alpha         : 1.-alpha     ) + phredConv.phred2Err[bq] / 4. );
    gls[7] *= ( phredConv.phred2Mat[bq] * (al == 0 ? alpha/2.      : 1.-alpha/2.  ) + phredConv.phred2Err[bq] / 4. );
    gls[8] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 0.0           : 1.0          ) + phredConv.phred2Err[bq] / 4. );    

    tmp = gls[0] + gls[1] + gls[2] + gls[3] + gls[4] + gls[5] + gls[6] + gls[7] + gls[8];
    gls[0] /= tmp;
    gls[1] /= tmp;
    gls[2] /= tmp;
    gls[3] /= tmp;
    gls[4] /= tmp;
    gls[5] /= tmp;
    gls[6] /= tmp;
    gls[7] /= tmp;
    gls[8] /= tmp;        
    //logdenom += log(tmp);
  }
  
  gls[0] += 1e-6;
  gls[1] += 1e-6;
  gls[2] += 1e-6;
  gls[3] += 1e-6;
  gls[4] += 1e-6;
  gls[5] += 1e-6;
  gls[6] += 1e-6;
  gls[7] += 1e-6;
  gls[8] += 1e-6;
  tmp = gls[0] + gls[1] + gls[2] + gls[3] + gls[4] + gls[5] + gls[6] + gls[7] + gls[8];  

  gls[0] /= tmp;
  gls[1] /= tmp;
  gls[2] /= tmp;
  gls[3] /= tmp;
  gls[4] /= tmp;
  gls[5] /= tmp;
  gls[6] /= tmp;
  gls[7] /= tmp;
  gls[8] /= tmp;          
  //logdenom += log(tmp);

  //return logdenom;
  return 0;  
}

double calculate_snp_droplet_GL(sc_snp_droplet_t* ssd, double* gls) {
  //double logdenom = 0;
  double tmp;
  gls[0] = gls[1] = gls[2] = 1.0;
  for(sc_snp_droplet_it_t it = ssd->begin(); it != ssd->end(); ++it) {
    uint8_t al = ( it->second >> 24 ) & 0x00ff;
    uint8_t bq = ( it->second >> 16 ) & 0x00ff;

    if ( al == 2 ) continue;

    gls[0] *= ((al==0) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0);
    gls[1] *= (0.5 - phredConv.phred2Err[bq]/3.0);
    gls[2] *= ((al==1) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0);
    tmp = gls[0] + gls[1] + gls[2];
    gls[0] /= tmp;
    gls[1] /= tmp;
    gls[2] /= tmp;
    //logdenom += log(tmp);
  }
  
  gls[0] += 1e-6;
  gls[1] += 1e-6;
  gls[2] += 1e-6;
  tmp = gls[0] + gls[1] + gls[2];
  gls[0] /= tmp;
  gls[1] /= tmp;
  gls[2] /= tmp;
  //logdenom += log(tmp);

  //return logdenom;
  return 0;
}
