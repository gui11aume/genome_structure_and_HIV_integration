public boolean isInteger( String input ) {
    try {
        Integer.parseInt( input );
        return true;
    }
    catch( NumberFormatException e ) {
        return false;
    }
}

public boolean isInteger( Integer input) {
   return true;
}

def print_help = {
    log.info ''
    log.info '   BHIVE Hi-C cooler pipeline  '
    log.info '-------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    cooler.nf [--resolution VALUE | --binfile BIN_FILE]  [optional arguments]'
    log.info ''
    log.info 'Computes Hi-C contacts.'
    log.info 'Datasets should be defined in "params.txt".'
    log.info ''
    log.info 'Required arguments:'
    log.info ''
    log.info '  either'
    log.info '    --resolution  Hi-C resolution (in bp).'
    log.info '  or'
    log.info '    --binfile     Hi-C bins in bed format.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '  --options     Path to alternative params file (default: params.txt)'
    log.info '  --out         Output path (default: .)'
    log.info ''
}

// Set defaults
params.help        = false
params.out         = '.'
params.options     = 'params.txt'
params.resolution  = null
params.binfile     = null
params.chromsizes  = null

if (params.help) {
   print_help()
   exit 0
}

if (!params.resolution && !params.binfile) {
   log.info 'error: required options --resolution or --binfile not set.'
   exit 1
}

if (params.resolution && params.binfile) {
   log.info 'error: set either --resolution or --binfile, but not both.'
   exit 1
}

binfile = Channel.create()
resolution = Channel.create()

if (!params.binfile) {
   if (!isInteger(params.resolution)) {
      log.info 'error: --resolution must be a numeric value.'
      exit 1
   } else {
      resolution << params.resolution
   }
} else {
   bins_path = file("${params.binfile}")
   if (!bins_path.isFile()) {
      log.info 'error: bin file not found'
      exit 1
   } else {
      binfile << [bins_path, bins_path.getBaseName()]
      binfile.close()
   }
}

resolution.close()

/**********************************************************************
**************************** PARSE  PARAMS ****************************
**********************************************************************/

/**********************************************************************
***************************** PARSE OPTIONS ***************************
**********************************************************************/
/*
** Mapping options format:
** 
** datasets:
** replicate,filename
** 
*/

def options = ['mapq': null, 'chromsize_path': null]

mfile = file("${params.options}")
if (!mfile.isFile()) {
   log.info "error: options file not found (${mfile})!"
   exit 1
}

// Varaibles
def status = 0

// Channels
file_getref = Channel.create()
datasets_   = Channel.create()

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
        fin = file("${t[-1]}")
        if (fin.isFile()) {
           datasets_ << [t[0],fin]
        } else {
           log.info "error: contact file not found, in '$line'."
           exit 1
        }
     } else {
        log.info "error: incorrect format of Hi-C contact dataset entry '$line' in $params.options"
        exit 1
     }
     break
   }
}

datasets_.close()

// Parse options
if (!options['mapq']) {
   log.info "error: mapq option not set in '${params.options}'."
   exit 1
}

if (!isInteger(options['mapq'])) {
   log.info "error: mapq value specified in '${params.options}' is not numeric!"
   exit 1
}

if (!options['chromsize_path']) {
   log.info "error: chromsize_path option not set in '${params.options}'."
   exit 1
}

csize_path = file("${options['chromsize_path']}")
if (!csize_path.isFile()) {
   log.info "error: chromsize file not found '${csize_path}'."
   exit 1
}

Channel.from(csize_path).into{ chromsizes_sort; chromsizes_bins }


// Group the same datasets.
datasets_.groupTuple().into{datasets}

/**********************************************************************
**************************** HI-C PIPELINE ****************************
**********************************************************************/

process coolerSort {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/sorted/", mode:'symlink'
   errorStrategy 'finish'

   // Cluster options
   cpus 8
   memory '24GB'

   input:
   set sample, file(reads) from datasets
   file chrfile from chromsizes_sort.first()

   output:
   set file("*.bgz"), file("*.bgz.tbi"), sample into sorted_contacts
   
   script:
   """
    cooler csort --sort-options '-S 2G --parallel=${task.cpus}' -i tabix -c1 1 -p1 2 -s1 3 -c2 4 -p2 5 -s2 6 -p ${task.cpus} -o hic_${sample}.scf.bgz <(zcat ${reads} | awk '{print \$2"\t"\$3"\t"\$4"\t"\$9"\t"\$10"\t"\$11}') ${chrfile}
   """
}

process coolerBins {
   // Process options
   publishDir path:"${params.out}/binfiles/", mode: 'symlink'
   errorStrategy 'finish'

   // Cluster options
   memory '8GB'
   time '1h'

   input:
   val binsize from resolution
   file chrfile from chromsizes_bins

   output:
   set file("*.bed"), binsize into binfile

   script:
   """
   cooler makebins -o bin_${binsize}.bed ${chrfile} ${binsize}
   """
}

process coolerLoad {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/coolerfiles/", mode: 'move'
   errorStrategy 'finish'
   
   // Cluster options
   cpus 8
   memory '24GB'
   time '8h'

   input:
   set file(bins), resolution from binfile.first()
   set file(contacts), file(tabix), sample from sorted_contacts

   output:
   file '*.cool' into cooler_files

   script:
   """
   cooler cload tabix -s 200 ${bins} ${contacts} ${sample}_${resolution}_q${options['mapq']}.cool
   """
}
