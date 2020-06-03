CREATE TABLE subject (
  subjid int(10) NOT NULL auto_increment,
  subjacc varchar(100),
  tissue varchar(100),
  tumortype varchar(100),
  index(subjid),
  index(subjacc),
  index(tumortype),
  index(tissue)
) ENGINE=MyISAM;

create table experiment (
  id int(10)  NOT NULL auto_increment,
  experiment varchar(100),
  dna_samplename varchar(100),
  rna_samplename varchar(100),
  index(id),
  index(experiment),
  index(dna_samplename),
  index(rna_samplename)
)  ENGINE=MyISAM;

create table clialab (
  driverstatus varchar(10),
  subjacc varchar(100),
  snptype varchar(10),
  chr varchar(10),
  pos int(20),
  snpid text,
  gene varchar(50),
  var varchar(200),
  altreadct int(10),
  readdepth int(10),
  maf float,
  callers varchar(200),
  exac_popmax_af float,
  impact varchar(20),
  numcosmicsubj int(10),
  index(driverstatus),
  index(subjacc),
  index(chr),
  index(pos), 
  index(gene),
  index(callers),
  index(impact)
) ENGINE=MyISAM;

CREATE TABLE samples (
  sid int(10) NOT NULL auto_increment,
  sample_id varchar(100),
  sample_name varchar(100),
  dnahybid int(10),
  rnahybid int(10),
  assay varchar(20),
  index(sid),
  index(sample_id),
  index(dnahybid),
  index(sample_name),
  index(assay)
) ENGINE=MyISAM;

CREATE TABLE seqgen (
   seqgen int(10) NOT NULL auto_increment,
   sid int(10),
   seqrunid varchar(200),
   demultiplexname varchar(200),
   lane int(10),
   clusters int(20),
   pfclusters int(20),
   percpf float,
   readyield double,
   percq30 float,
   avgqualscore float,
   perclane float,
   dempstatus enum('pass','fail'),
   index(seqgen),
   index(sid),
   index(seqrunid),
   index(lane),
   index(dempstatus),
   index(samplename)
) ENGINE=MyISAM;


CREATE TABLE demultiplex (
   sid int(10),
   flowcellid varchar(200),
   samplename varchar(200),
   lane int(10),
   clusters int(20),
   pfclusters int(20),
   percpf float,
   readyield double,
   percq30 float,
   avgqualscore float,
   perclane float,
   dempstatus enum('pass','fail'),
   index(sid),
   index(flowcellid),
   index(lane),
   index(dempstatus),
   index(samplename)
) ENGINE=MyISAM;

CREATE TABLE alignment (
  sid int(10),      
  samplename varchar(200),
  totalreads int(20),
  maprate float,
  propair float,
  percontarget float,
  meanmapq float,
  percqualread float,
  mediantmismatch float,
  indelrate float,
  errorrate float,
  uniqreads int(10),
  percdups float,
  medianinsert int(10),
  meaninsert int(10),
  stdinsert int(10),
  avg_depth int(10),
  median_depth int(10),
  perc100X float,
  alignstatus enum('pass','fail'),
  alignmentdate date,
  workflow varchar(10),
  initiator varchar(20),
  index(sid),
  index(samplename),
  index(alignstatus),
  index(workflow),
  index(alignmentdate),
  index(initiator)
) ENGINE=MyISAM;  

create table seqrun (
  seqrunid varchar(100),
  flowcellid varchar(200),
  seqmachine varchar(50),
  phixerror float,
  phixswitchrate float,
  rundate date,
  runtech varchar(100),
  index(flowcellid),
  index(seqmachine),
  index(runtech)
) ENGINE=MyISAM;

create table coverage (
  sid int(10),      
  samplename varchar(200),
  chrom varchar(10),
  cstart int(10),
  cend int(10),
  exonname text,
  mindepth int(10),
  maxdepth int(10),
  mediandepth int(10),
  avgdepth int(10),
  frac100plus int(10),
  bp100plus int(10),
  totalbp int(10),
  tier1exon int(1),
  index(sid),
  index(samplename)  
) ENGINE=MyISAM;  

CREATE TABLE nucextract (
  nucextractid int(10) NOT NULL auto_increment,
  subjacc varchar(100),
  extractname varchar(100),
  tclass varchar(100),
  fixtype varchar(10),
  specimen varchar(100),
  thickness float,
  numslides int(10),
  extracttech varchar(100),
  dateextract date,
  index(nucextractid),
  index(subjacc),
  index(tclass),
  index(fixtype),
  index(specimen),
  index(extracttech)
) ENGINE=MyISAM;

CREATE TABLE dnaextract (
  dnaextractid int(10) NOT NULL auto_increment,
  nucextractid int(10),
  extractdin float,
  dnayield float,
  extractstatus enum('pass','fail'),
  index(dnaextractid),
  index(nucextractid),
  index(extractstatus)
) ENGINE=MyISAM;

CREATE TABLE rnaextract (
  rnaextractid int(10) NOT NULL auto_increment,
  nucextractid int(10),
  perc200ntplus float,
  rnayield float,
  extractstatus enum('pass','fail'),
  index(rnaextractid),
  index(nucextractid),
  index(extractstatus)
) ENGINE=MyISAM;

create table dnalibprep (
  dnalibprepid int(10) NOT NULL auto_increment,
  dnaextractid int(10),
  libraryname varchar(50),
  inputdna float,
  meanlibfragsize int(10),
  libyield float,
  perc200ntplus float,
  libstatus enum('pass','fail'),
  dnahybid int(10),
  bcindex varchar(20),
  libtech varchar(100),
  datelib date,
  index(dnalibprepid),
  index(libraryname),
  index(dnaextractid),
  index(dnahybid),
  index(bcindex),
  index(libstatus),
  index(libtech),
  index(datelib)
) ENGINE=MyISAM;

create table rnalibprep (
  rnalibprepid int(10) NOT NULL auto_increment,
  rnaextractid int(10),
  libraryname varchar(50),
  perc200ntplus float,
  inputlib float,
  libstatus enum('pass','fail'),
  prehybsize int(10),
  inputhyb float,
  posthybsize int(10),
  hybconcentration float,
  hybstatus enum('pass','fail'),
  libtech varchar(100),
  datelib date,
  index(rnalibprepid),
  index(libraryname),
  index(rnaextractid),
  index(libstatus),
  index(hybstatus),
  index(libtech),
  index(datelib)
) ENGINE=MyISAM;

create table dnahyb (
  dnahybid int(10) NOT NULL auto_increment,
  hybmolarity float,
  meanfragsize int(10),
  hybyield float,
  perc200ntplus float,
  hybstatus enum('pass','fail'),
  hybtech varchar(100),
  datehyb date,
  index(dnahybid),
  index(hybstatus),
  index(hybtech),
  index(datehyb)
) ENGINE=MyISAM;

CREATE TABLE rnaqm (
  sid int(10),
  rin float,
  extractyield float,
  inputrna float,
  perc200ntplus float,
  prelibconcentration float,
  avgprehyblibfragsize float,
  hybinput int(10),
  avgposthyblibfragsize int(20),
  posthybyield float,
  index(sid)
) ENGINE=MyISAM;

