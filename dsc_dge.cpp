#include "dsc_dge.h"
#include "tsv_reader.h"

dsc_dge::dsc_dge(const char* bcdf, const char* genef, const char* mtxf) {
  notice("dsc_dge(%s, %s, %s)", bcdf, genef, mtxf);
  notice("Reading %s",bcdf);
  tsv_reader tsv_bcdf(bcdf);
  for(int32_t i=0; tsv_bcdf.read_line() > 0; ++i) {
    bcds.push_back(tsv_bcdf.str_field_at(0));
    bcd2idx[bcds.back()] = i;
  }
  tsv_bcdf.close();

  notice("Reading %s",genef);  
  tsv_reader tsv_genef(genef);
  for(int32_t i=0; tsv_genef.read_line() > 0; ++i) {
    gene_ids.push_back(tsv_genef.str_field_at(0));
    if ( tsv_genef.nfields == 1 ) {
      gene_symbs.push_back(tsv_genef.str_field_at(0));      
    }
    else {
      gene_symbs.push_back(tsv_genef.str_field_at(1));      
    }
  }
  tsv_genef.close();

  notice("Reading %s",mtxf);    
  tsv_reader tsv_mtxf(mtxf);
  int64_t total_umis = 0, sum_umis = 0;
  bcd_umis.resize(bcds.size(), 0);
  gene_umis.resize(gene_ids.size(), 0);
  for(int32_t i=0; tsv_mtxf.read_line() > 0; ++i) {
    if ( i % 1000000 == 0 ) notice("Reading %d lines in %s",i,mtxf);
    if ( i < 2 ) {
      if ( tsv_mtxf.str_field_at(0)[0] != '%' )
	error("First two lines are expected to start with %%, but observed %s", tsv_mtxf.str_field_at(0));
    }
    else if ( i == 2 ) {
      if ( tsv_mtxf.int_field_at(0) != (int32_t)gene_ids.size() )
	error("Inconsistency in DGE matrix %d genes are observed in %s, but %d are recorded in %s", (int32_t)gene_ids.size(), genef, tsv_genef.int_field_at(0), mtxf);      
      if ( tsv_mtxf.int_field_at(1) != (int32_t)bcds.size() )
	error("Inconsistency in DGE matrix %d barcodes are observed in %s, but %d are recorded in %s", (int32_t)bcds.size(), bcdf, tsv_mtxf.int_field_at(1), mtxf);
      total_umis = tsv_mtxf.int_field_at(2);
    }
    else {
      if ( tsv_mtxf.nfields != 3 ) {
	error("Expecting 3 columns but observed %d at line %d", tsv_mtxf.nfields, i);
      }
      
      int32_t igene = tsv_mtxf.int_field_at(0)-1;
      int32_t ibcd  = tsv_mtxf.int_field_at(1)-1;
      int32_t umi   = tsv_mtxf.int_field_at(2);
      bcd_gene_umi[ibcd][igene] = umi;
      gene_bcd_umi[igene][ibcd] = umi;
      bcd_umis[ibcd]  += umi;
      gene_umis[igene] += umi;
      sum_umis += umi;
    }
  }
  tsv_mtxf.close();

  notice("Successfully loaded %lld UMIs across %d genes and %d barcodes from %s", total_umis, (int32_t)gene_ids.size(), (int32_t)bcds.size(), mtxf);
  if ( total_umis != sum_umis )
    warning("Loaded %lld UMIs across %d genes and %d barcodes from %s, but it is different from expected value %lld", sum_umis, (int32_t)gene_ids.size(), (int32_t)bcds.size(), mtxf, total_umis);
}
