function genome_info=init_genome(config_fname)
% genome_info=init_genome(config_fname)  
%
% global g_genome_info ;

global g_genome_info genome_strings;

genome_strings = [];
g_est_fnames={} ;
g_prot_fnames={} ;
g_cdna_fnames={} ;
g_flcdna_fnames={} ;
g_contig_names={} ;
g_flat_fnames={} ;
g_fasta_fnames={} ;
g_annotation_fnames={} ;
g_alphabet='acgt' ;
g_basedir='' ;
g_max_intron_len=1000000 ;
g_min_intron_len=10 ;
g_merge_est_transcripts=1 ;

[fd msg]=fopen(config_fname, 'r') ;
if fd<1, error('could not open file %s: %s',config_fname,msg); end
while ~feof(fd),
  line = fgetl(fd) ;
  if ~ischar(line), break ; end ;
  if isempty(line), continue ; end ;
  if line(1)=='#', continue ; end ;
  if length(line)>8 && isequal(line(1:7),'BASEDIR')
    g_basedir=deblank(line(9:end)) ;
  elseif length(line)>21 && isequal(line(1:21), 'MERGE_EST_TRANSCRIPTS')
    g_merge_est_transcripts = str2num(line(22:end)) ; 
  elseif length(line)>15 && isequal(line(1:15),'MAX_DOUBLE_QGAP')
  elseif length(line)>15 && isequal(line(1:15),'MAX_DOUBLE_TGAP')
  elseif length(line)>17 && isequal(line(1:17),'MAX_CDNA_DELETION')
  elseif length(line)>18 && isequal(line(1:18),'MAX_CDNA_INSERTION')
  elseif length(line)>21 && isequal(line(1:21),'TERMINAL_EXON_END_TOL')
  elseif length(line)>14 && isequal(line(1:14),'MIN_INTRON_LEN')
    g_min_intron_len = str2num(line(15:end)) ; 
    assert(~isempty(g_min_intron_len) & g_min_intron_len<1000) ;
  elseif length(line)>14 && isequal(line(1:14),'MAX_INTRON_LEN')
    g_max_intron_len = str2num(line(15:end)) ; 
    assert(~isempty(g_max_intron_len) & g_max_intron_len>1000) ;
  elseif length(line)>18 && isequal(line(1:18),'MIN_EST_COVER_FRAC')

  elseif length(line)>19 && isequal(line(1:19),'MIN_PROT_COVER_FRAC')

  elseif length(line)>19 && isequal(line(1:19),'MIN_CDNA_COVER_FRAC')
  
  elseif length(line)>18 && isequal(line(1:18),'BLAT_BEST_HIT_ONLY')
  
  elseif length(line)>20 && isequal(line(1:20),'BLAT_BEST_HIT_MARGIN')

  elseif length(line)>8 && isequal(line(1:8),'ALPHABET')
    g_alphabet=deblank(line(10:end)) ;
  elseif length(line)>7 && isequal(line(1:7),'CONTIGS'),
    for i=1:str2num(line(8:end))
      l2=fgetl(fd) ; 
      idx=strfind(l2,'  ') ;
      while ~isempty(idx), 
        l2(idx)=[] ;
        idx=strfind(l2,'  ') ;
      end ;
      l2(l2==' ')=sprintf('\t') ;
      elems=separate(l2) ; 
      assert(length(elems)==3) ;
      g_contig_names{i}=elems{1} ;
      g_flat_fnames{i}=path_with_basedir(elems{2}, g_basedir) ;
      g_fasta_fnames{i}=path_with_basedir(elems{3}, g_basedir) ; 
    end ;
  elseif length(line)>8 && isequal(line(1:8),'ESTFILES'),
    for i=1:str2num(line(9:end))
      g_est_fnames{i}=path_with_basedir(deblank(fgetl(fd)), g_basedir) ; 
    end ;
  elseif length(line)>9 && isequal(line(1:9),'PROTFILES'),
    for i=1:str2num(line(10:end))
      g_prot_fnames{i}=path_with_basedir(deblank(fgetl(fd)), g_basedir) ; 
    end ;
  elseif length(line)>9 && isequal(line(1:9),'CDNAFILES'),
    for i=1:str2num(line(10:end))
      g_cdna_fnames{end+1}=path_with_basedir(deblank(fgetl(fd)), g_basedir) ; 
    end ;
  elseif length(line)>11 && isequal(line(1:11),'FLCDNAFILES'),
    for i=1:str2num(line(12:end))
      g_flcdna_fnames{end+1}=path_with_basedir(deblank(fgetl(fd)), g_basedir) ; 
    end ;
  elseif length(line)>15 && isequal(line(1:15),'ANNOTATIONFILES'),
    for i=1:str2num(line(16:end))
      g_annotation_fnames{i}=path_with_basedir(deblank(fgetl(fd)), g_basedir) ; 
    end ;
  elseif ~isempty(line),
    error(sprintf('cannot understand "%s" in config file', line)) ;
  end ;
end ;
fclose(fd) ;

genome_info=[] ;
genome_info.est_fnames=g_est_fnames ;
genome_info.prot_fnames=g_prot_fnames ;
genome_info.cdna_fnames=g_cdna_fnames ;
genome_info.flcdna_fnames=g_flcdna_fnames ;
genome_info.contig_names=g_contig_names ;
genome_info.flat_fnames=g_flat_fnames ;
genome_info.fasta_fnames=g_fasta_fnames ;
genome_info.annotation_fnames=g_annotation_fnames ;
genome_info.alphabet=g_alphabet ;
genome_info.basedir=g_basedir ;
genome_info.max_intron_len = g_max_intron_len ;
genome_info.min_intron_len = g_min_intron_len ;
genome_info.merge_est_transcripts = g_merge_est_transcripts ;
g_genome_info=genome_info ;

return ;


function p=path_with_basedir(p1,basedir)
  
  if p1(1)=='/',
    p=p1; 
    return ;
  else
    if basedir(end)=='/',
      p=[basedir p1] ;
    else
      p=[basedir '/' p1] ;
    end ;
  end ;
  
  
