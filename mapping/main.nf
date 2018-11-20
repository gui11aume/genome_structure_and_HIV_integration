public boolean isInteger( String input ) {
    try {
        Integer.parseInt( input );
        return true;
    }
    catch( NumberFormatException e ) {
        return false;
    }
}

def print_help = {
    log.info ''
    log.info '  BHIVE Hi-C mapping pipeline  '
    log.info '-------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    hic-mapping.nf --index BWA_INDEX [optional arguments]'
    log.info ''
    log.info 'Computes Hi-C contacts.'
    log.info 'Datasets should be defined in "params.txt".'
    log.info ''
    log.info 'Required arguments:'
    log.info ''
    log.info '  --index       Path to BWA index file (required)'
    log.info ''
    log.info 'Optional arguments:'
    log.info '  --options     Path to alternative mapping params file (default: params.txt)'
    log.info '  --out         Output path (default: .)'
    log.info ''
}

src_dir = "src"
bin_dir = "bin"

hic_parser = Channel.fromPath("${bin_dir}/parse_contacts")

// Set defaults
params.help        = false
params.out         = '.'
params.options     = 'params.txt'
params.index       = null

if (params.help) {
   print_help()
   exit 0
}

if (!file("${bin_dir}/parse_contacts").isFile()) {
   log.info "parser file, expected in '${bin_dir}/parse_contacts' not found."
   exit 1
}

/**********************************************************************
**************************** PARSE  PARAMS ****************************
**********************************************************************/


bwa_index   = Channel.create()
index_files = Channel.create()
// 0. Find BWA index
if (params.index) {
   // Check bwa index files.
   fasta_ref = file("${params.index}.bwt")
   index_ref = file("${params.index}")
   if (fasta_ref.isFile()) {
      if (!index_ref.isFile()) {
         index_ref.toFile().createNewFile()
      }
      bwa_index << file("${params.index}")
      index_files << file("${index_ref}.{bwt,amb,ann,pac,sa}")
   } else if (index_ref.getExtension() in ['bwt','amb','ann','pac','sa'] && index_ref.isFile()) {
      index_str = (index_ref.getParent().equals(null) ? './' : index_ref.getParent()) + index_ref.getBaseName()
      index_ref = file(index_str)
      if (!index_ref.isFile()) {
         index_ref.toFile().createNewFile()
      }
      bwa_index << index_ref
      index_files << file("${index_ref}.{bwt,amb,ann,pac,sa}")
   } else {
      log.info "error: BWA index not found in '${params.index}'."
      exit 1
   }
} else {
   log.info "error: '--index' option not specified."
   exit 1
}
bwa_index.close()
index_files.close()


/**********************************************************************
***************************** PARSE OPTIONS ***************************
**********************************************************************/
/*
** Mapping options format:
** 
** datasets:
** replicate,filename/URL/{SRR,ERR,DRR} reference
** 
*/

def options = ['re_db_path': null, 're_organism': null, 're_name': null, 'mapq': null, 'ins_size': null]

mfile = file("${params.options}")
if (!mfile.isFile()) {
   log.info "error: options file not found (${mfile})!"
   exit 1
}

// Varaibles
def status = 0

// Channels
file_getref = Channel.create()
datasets    = Channel.create()

// Parse file by lines.
mfile.eachLine { line ->
   if (line =~ /^#/ || line =~ /^\s*$/) {
      return null
   } else if (line =~ /^datasets:/) {
      status = 1
      return null
   }
   switch (status) {
   case 0:
     t = line.split(/\s*=\s*/).toList()
     if (!options.containsKey(t[0])) {
        log.info "unknown option: '${t[0]}' in $params.options"
        exit 1
     }
     options[t[0]] = t[1]
     break
   case 1:
     t = line.split(/\s*,\s*/).toList()
     if (t.size() == 2) {
        if (t[1] =~ /^SRR|^ERR|^DRR/) {
           file_getref << [[t[-1]],t[0]]
        } else {
           log.info "error: dataset iPCR entries must specify 2 files, 2 URL or a GEO reference starting with SRR/ERR/DRR. Error in entry '$line'"
           exit 1
        }
     } else if (t.size() == 3) {
        if (t[1] =~ /^http|^ftp/ && t[2] =~ /^http|^ftp/) {
           file_getref << [[t[-2],t[-1]],t[0]]
        } else {
           read1 = file("${t[-2]}")
           read2 = file("${t[-1]}")
           if (read1.isFile() && read2.isFile()) {
              datasets << [[read1,read2],t[0]]
           } else {
              log.info "error: iPCR files not found, in '$line'. (URLs must start with http/ftp)"
              exit 1
           }
        }
     } else {
        log.info "error: incorrect format of iPCR dataset entry '$line' in $params.options"
        exit 1
     }
     break
   }
}

// Parse options
if (!options['re_db_path'] || !options['re_organism'] || !options['re_name'] || !options['mapq'] || !options['ins_size']) {
   log.info "error: 're_db_path', 're_organism', 're_name', 'mapq' and 'ins_size' options must be defined in $params.options before 'datasets:'."
   exit 1
}

if (!isInteger(options['mapq']) || !isInteger(options['ins_size'])) {
   log.info "error: 'mapq' and 'ins_size' options defined in $params.options must be numeric."
   exit 1
}

// DB path
db_path = Channel.fromPath(options['re_db_path'])

// Group the same references/URLS to avoid redundant downloads.
file_getref.close()
file_getref.groupTuple().into{gfile_ref}

process getDataset {
   // Process options
   tag "${data[0]}_${data[1]}"
   publishDir path:"${params.out}/datasets/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val data from gfile_ref
   output:
   set file('*.fastq.gz'), replicate into datasets
   script:
   replicate = data[1]
   ref = data[0]
   if (ref.size() == 1) {
      if (ref[0] =~ /^SRR|^ERR|^DRR/) {
         """
         fastq-dump --split-files --gzip -A ${ref[0]}
         rm -f ~/ncbi/public/sra/${ref[0]}.sra
         """
      } else {
         log.info "error: only one Hi-C read specified (must be PE files or sigle GEO reference)!"
      }
   } else if (ref.size() == 2) {
      """
      wget ${ref[0]}
      wget ${ref[1]}
      """
   }
}

/**********************************************************************
**************************** HI-C PIPELINE ****************************
**********************************************************************/

process mapContacts {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/mapped/", mode:'symlink'
   errorStrategy 'finish'

   // Cluster options
   cpus 12
   memory '32GB'
   time '7d'

   input:
   set file(reads), sample from datasets
   file index_path from bwa_index.first()
   file index_files from index_files.first()

   output:
   set file("*.bam"), sample into raw_contacts
   
   script:
   """
    OPTS_36='-B3 -O5 -P -k17 -U0 -L0,0 -T23'
    OPTS_40='-P -k17 -U0 -L0,0 -T25'
    OPTS_50='-P -k17 -U0 -L0,0 -T25'
    SEQLEN=\$((\$(zcat -f ${reads[0]} | head -2 | tail -1 | wc -c) - 1));
    if [ \$SEQLEN -le 30 ]; then \
      OPTS=\$OPTS_26; \
    elif [ \$SEQLEN -le 39 ]; then \
      OPTS=\$OPTS_36; \
    elif [ \$SEQLEN -le 46 ]; then \
      OPTS=\$OPTS_40; \
    else \
      OPTS=\$OPTS_50; \
    fi;
    echo \$SEQLEN
    echo \$OPTS
    echo ${index_path}
    bwa mem -t ${task.cpus} \$OPTS ${index_path} ${reads[0]} ${reads[1]} | samtools view -bS - > hic_${sample}.bam
   """
}

process computeContacts {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/contacts/", mode:'move'

   // Cluster options
   cpus 12
   memory '16GB'

   input:
   set file(bam), sample from raw_contacts
   file parser from hic_parser.first()
   file redb from db_path.first()

   output:
   set file("*.hcf"), sample into contact_files
   
   script:
   """
   chmod +x ${parser}
   ./${parser} ${options['re_organism']} ${options['re_name']} <(samtools view ${bam}) ${options['mapq']} ${options['ins_size']} > hic_${sample}.hcf
   """
}
