# The GATK Best Practices Workflow from Github (Oct 16, 2020 version 1.0)
# v1.1b
# Usage: source the script
#

# TODO: 1, check input file existance
#       5, append gatk Funcotator to variants calling function

#       7, add germline CNV calling process

# each process must have a standard name and header (sample.name.str, temp.str, output, sleep.time etc.)

library(parallel)
library(tools)
library(stringr)

options(warn = -1)  # suppress all warnings

Process.create <- function(file, samples, steps) {
   process.df = data.frame(matrix('', nrow = length(samples), ncol = length(steps)))
   rownames(process.df) = samples
   colnames(process.df) = steps
   write.table(process.df, file, row.names = T, col.names = T, sep = '\t', quote = F, na = '')
}

Process.status.change <- function(file, sample, step, change.to = 'Done') {
   process.df = read.csv(file, row.names = 1, sep = '\t', stringsAsFactors = F)
   process.df[sample, step] = change.to
   write.table(process.df, file, row.names = T, col.names = T, sep = '\t', quote = F, na = '')
}

Process.status.add <- function(file, sample, step, add = 1) {
   process.df = read.csv(file, row.names = 1, sep = '\t', stringsAsFactors = F)
   content.str = process.df[sample, step]
   number.c = str_split(content.str, '/', simplify = T)[1,]
   number.c[1] = as.integer(number.c[1]) + 1
   process.df[sample, step] = paste(number.c, collapse = '/')
   write.table(process.df, file, row.names = T, col.names = T, sep = '\t', quote = F, na = '')
}

get_machine_mem <- function(floor = T) {
   library(stringr)
   temp = file('/proc/meminfo', open = 'rt')
   mem.str = readLines(temp)[1]
   mem.int = as.integer(str_extract(mem.str, pattern = '\\d+')) / 1024 / 1024
   close(temp)
   return(floor(mem.int))
}
get_cpu_count <- function() {
   library(stringr)
   temp = file('/proc/cpuinfo', open = 'rt')
   cpu.c = readLines(temp)
   close(temp)
   cpu.c = cpu.c[str_starts(cpu.c, pattern = 'processor')]
   return(length(cpu.c))
}
get_fastq_read_length <- function(fastq.file) {
   library(stringr)
   if (str_detect(fastq.file, pattern = '.gz$')) {
      fq = gzfile(fastq.file, open = 'rt')
   } else {
      fq = file(fastq.file, open = 'rt')
   }
   temp = readLines(con = fq, n = 4, ok = F, warn = F)
   close(fq)
   return(str_length(temp[2]))

}

# get base name such as 1801169
get_sample_name <- function(file) {
   library(tools)
   temp.str = file_path_sans_ext(basename(file))
   sample.name.str = str_extract(temp.str, pattern = '\\S+(?=-N|-T|-DN|-DT|-RN|-RT)')
   return(sample.name.str)
}


# fastq must be '\\S+(-T|-DT|-RT|-N|-DN|-RN)_R(1|2)\\.fq\\.gz$'
generate.fastq.list <- function(fastq.list.file) {
   library(stringr)
   library(tools)

   # read fastq.list
   fastq.files.c = read.csv(fastq.list.file, sep = '\t', header = F, stringsAsFactors = F, comment.char = '#', strip.white = T, blank.lines.skip = T)$V1

   # 994278-DT_R2.fq.gz
   file.pattern = '\\S+(-T|-DT|-RT|-N|-DN|-RN)_R(1|2)\\.fq\\.gz$'
   sample.name.pattern = '\\S+(?=_R(1|2))'

   fastq.list = list()
   for (i in 1:length(fastq.files.c)) {
      if (!str_detect(fastq.files.c[i], pattern = file.pattern)) {
         message(fastq.files.c[i], 'Unknown file. Skip')
         next
      }  # skip the non fastq file

      temp.str = file_path_sans_ext(basename(fastq.files.c[i]))  # 1565230-DT_R2.fq
      sample.name.str = str_extract(temp.str, pattern = sample.name.pattern)  # 1565230-DT
      if (str_detect(temp.str, '_R1')) {
         type = 'R1'
      } else if (str_detect(temp.str, '_R2')) {
         type = 'R2'
      } else {
         message(fastq.files.c[i], ' Unknown type. Skip')
         next
      }

      if (sample.name.str %in% names(fastq.list)) {
         fastq.list[[sample.name.str]][type] = fastq.files.c[i]
      } else {
         fastq.list[[sample.name.str]] = structure(c(fastq.files.c[i]), names = type)
      }
   }
   return(fastq.list)
}


# bam file must be '\S+(-T|-DT)\.bam$'
# bams.list.file is a vector cotaining all bam files or a text file in which each line contains a bam file
# allow.type is one of 'both', 'tumor' and 'normal'
# return.type is either 'list' or 'vector'
generate.bam.list <- function(bams.list, allow.type = 'both', return.type = 'list') {
   library(stringr)
   library(tools)

   # read bam files
   if (length(bams.list) == 1) {
      bam.files.c = normalizePath(read.csv(bams.list, sep = '\t', header = F, stringsAsFactors = F, comment.char = '#', strip.white = T, blank.lines.skip = T)$V1, mustWork = T)
   } else {
      bam.files.c = bams.list
   }

   if (allow.type == 'both') {
      file.pattern = '\\S+(-T|-DT|-RT|-N|-DN|-RN)\\.bam$'  # 1801167-DT.bam
      sample.name.pattern = '\\S+(?=-T|-DT|-RT|-N|-DN|-RN)'  # 1801167

   } else if (allow.type == 'normal') {
      file.pattern = '\\S+(-N|-DN|-RN)\\.bam$'
      sample.name.pattern = '\\S+(?=-N|-DN|-RN)'

   } else if (allow.type == 'tumor') {
      file.pattern = '\\S+(-T|-DT|-RT)\\.bam$'
      sample.name.pattern = '\\S+(?=-T|-DT|-RT)'
   }

   type.tumor.pattern = '(-T|-DT|-RT)'
   type.normal.pattern = '(-N|-DN|-RN)'

   filterd.files.c = c()
   for (file in bam.files.c) {
      if (str_detect(file, pattern = file.pattern)) {
         filterd.files.c = c(filterd.files.c, file)
      }
   }

   if (return.type == 'vector') {
      return(filterd.files.c)
   }

   bam.list = list()
   for (i in 1:length(filterd.files.c)) {

      temp.str = tools::file_path_sans_ext(basename(filterd.files.c[i]))  # 1801167-DT
      sample.name.str = str_extract(temp.str, pattern = sample.name.pattern)  # 1801167
      if (str_detect(temp.str, pattern = type.tumor.pattern)) {
         type = 'tumor.bam'
      } else if (str_detect(temp.str, pattern = type.normal.pattern)) {
         type = 'normal.bam'
      } else {
         message('Unknown type. Skip')
         next
      }

      if (sample.name.str %in% names(bam.list)) {
         bam.list[[sample.name.str]][type] = filterd.files.c[i]
      } else {
         bam.list[[sample.name.str]] = structure(c(filterd.files.c[i]), names = type)
      }
   }
   return(bam.list)

}


# command.c = c(prefer = '', fail.safe = '')
run <- function(command.c, message, quit.on.fail = T) {
   Sys.sleep(runif(1, 0.1, 2))
   message(sprintf('\n[%s] %s', Sys.time(), message))
   r = system(command.c['prefer'], ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.c['prefer']))
      return(0)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.c['prefer']))
      message(sprintf('\nTry fail-safe: %s', command.c['fail.safe']))
      r = system(command.c['fail.safe'], ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.c['fail.safe']))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.c['fail.safe']))
         if (quit.on.fail) {
            quit('no', status = r)
         } else {
            return(1)
         }
      }
   }
}

run.single.command <- function(command, sample.name, abbr, verbose = T) {
   if (verbose) {message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name, abbr))}
   r = system(command, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      if (verbose) {message(sprintf('\n[%s] %s: %s (SUCCESS)', Sys.time(), sample.name, abbr))}
   } else {
      message(sprintf('\n[%s] %s: %s (Fail %s: %s)', Sys.time(), sample.name, abbr, r, command))
   }
   return(r)
}

FAT.NODES = c(s001 = '192.168.11.31',
              s002 = '192.168.11.32')

THIN.NODES = c(c001 = '192.168.11.11',
               c002 = '192.168.11.12',
               c003 = '192.168.11.13',
               c004 = '192.168.11.14',
               c005 = '192.168.11.15',
               c006 = '192.168.11.16',
               c007 = '192.168.11.17',
               c008 = '192.168.11.18')

ssh.apply <- function(commands, nodes) {

}





# run in parallele mode
# temp.folder = sprintf('temp/preprocess.dna/%s', sample.name.str)
# named.vector = c(R1 = '~/projects/test/fastq/1801169-DT_R1.fq', R2 = '~/projects/test/fastq/1801169-DT_R2.fq')
# bwa = 'bwa'
# samtools = 'samtools'
# gatk = '~/bin/gatk-4.1.8.1/gatk'
# scatter.count = 8
# ref.fasta = '~/sda1/human_g1k_v37/human_g1k_v37.fasta'
# ref.dict = '~/sda1/human_g1k_v37/human_g1k_v37.dict'
# known.sites.VCFs = c('~/support_file/b37/dbsnp_138.b37.vcf', '~/support_file/b37/Mills_and_1000G_gold_standard.indels.b37.vcf')
# output: {output.folder}/{sample.name.str}.bam
..Preprocessing.DNA <- function(named.vector, mapper, mapper.index, samtools, gatk, scatter.count, keep.unBQSR, ref.fasta, ref.dict, known.sites.VCFs) {
   Sys.sleep(runif(1, 0.1, 2))

   fastq.1.str = named.vector['R1']
   fastq.2.str = named.vector['R2']
   mapper = mapper
   mapper.index = mapper.index
   samtools = samtools
   gatk = gatk
   scatter.count = scatter.count
   ref.fasta = ref.fasta
   ref.dict = ref.dict
   known.sites.VCFs = known.sites.VCFs

   sample.name.str = str_extract(basename(fastq.1.str), pattern = '^\\S+(?=_R1)')  # 1801169-PT
   if (sample.name.str != str_extract(basename(fastq.2.str), pattern = '^\\S+(?=_R2)')) {
      message(sprintf('%s and %s have different sample names.', fastq.1.str, fastq.2.str))
      return(1)
   }

   temp.folder = sprintf('temp/preprocess.dna/%s', sample.name.str)  # temp/preprocess.dna/1801169-PT
   unlink(temp.folder, force = T, recursive = T)
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = 'bamfiles'
   dir.create(output.folder, showWarnings = F)

   CreateSequenceGroupingTSV <- function(ref.dict) {
      library(stringr)
      result.list = list()

      dict.df = read.csv(ref.dict, header = F, sep = '\t', stringsAsFactors = F)
      dict.df = dict.df[dict.df$V1 == '@SQ', ]
      dict.df$chr = str_sub(dict.df$V2, start = 4)
      dict.df$len = as.integer(str_sub(dict.df$V3, start = 4))
      dict.df = dict.df[order(dict.df$len, decreasing = T), ]
      rownames(dict.df) = NULL
      longest.length.int = max(dict.df$len)  # 249250621 for hg19

      group.c = c()
      current.chrom.str = dict.df$chr[1]
      current.length.int = dict.df$len[1]
      for (i in 2:nrow(dict.df)) {
         if (current.length.int + dict.df$len[i] >= longest.length.int) {
            group.c = c(group.c, current.chrom.str)
            current.chrom.str = dict.df$chr[i]
            current.length.int = dict.df$len[i]
         } else {
            current.chrom.str = paste0(current.chrom.str, ' -L ', dict.df$chr[i])
            current.length.int = current.length.int + dict.df$len[i]
         }
      }
      group.c = c(group.c, current.chrom.str)
      group.c = paste('-L', group.c)

      result.list[['mapped']] = group.c
      result.list[['unmapped']] = c(group.c, '-L unmapped')

      return(result.list)

   }

   # ==========================map               output:{temp.folder}/querysort.bam==========================
   if (str_detect(mapper, 'bwa')) {
      command.str = sprintf("%s mem -K 100000000 -P -t %s -o %s/raw.sam -R '@RG\\tID:%s\\tSM:%s\\tPL:Illumina' %s %s %s", mapper, scatter.count, temp.folder, sample.name.str, sample.name.str, mapper.index, fastq.1.str, fastq.2.str)
   } else if (str_detect(mapper, 'minimap2')) {
      command.str = sprintf("%s -t %s -a -o %s/raw.sam -R '@RG\\tID:%s\\tSM:%s\\tPL:Illumina' -x sr -L --MD %s %s %s", mapper, scatter.count, temp.folder, sample.name.str, sample.name.str, mapper.index, fastq.1.str, fastq.2.str)
   } else {
      message(sprintf('Mapper must be either bwa mem or minimap2 instead of %s', mapper))
      return(1)
   }

   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'map'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = F)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/raw.sam

   # command.str = sprintf('%s view -1 -bh -o %s/unsorted.bam %s/raw.sam', samtools, temp.folder, temp.folder)
   # message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'convert to BAM'))
   # r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   # if (r == 0) {
   #    message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   # } else {
   #    message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
   #    quit('no', status = r)
   # }
   # # output {temp.folder}/unsorted.bam

   #command.str = sprintf('%s SortSam -I %s/unsorted.bam -O %s/querysort.bam -SO queryname', gatk, temp.folder, temp.folder)
   command.str = sprintf('%s sort -@%s -n -o %s/querysort.bam %s/raw.sam', samtools, scatter.count, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'sort by queryname'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }

   # output {temp.folder}/querysort.bam

   # ==========================MarkDuplicates    output:{temp.folder}/dup_rem.bam==========================
   # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
   # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
   # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
   command.str = sprintf("%s --java-options -Xms8G MarkDuplicates -I %s/querysort.bam -O %s/dup_rem.bam --METRICS_FILE %s/MarkDuplicates.metrics.txt --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER 'queryname'", gatk, temp.folder, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MarkDuplicates'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output:
   # 1, {temp.folder}/dup_rem.bam
   # 2, {temp.folder}/MarkDuplicates.metrics.txt

   # ==========================SortAndFixTags    output:{temp.folder}/fix.bam==========================
   command.str = sprintf("%s --java-options -Xms8G SortSam -I %s/dup_rem.bam -O %s/coordinatesort.bam --SORT_ORDER 'coordinate' --CREATE_INDEX false --CREATE_MD5_FILE false", gatk, temp.folder, temp.folder )
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'sort by coordinate'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output:  {temp.folder}/coordinatesort.bam

   command.str = sprintf('%s --java-options -Xms8G SetNmMdAndUqTags -I %s/coordinatesort.bam -O %s/fix.bam --CREATE_INDEX true --CREATE_MD5_FILE true --REFERENCE_SEQUENCE %s', gatk, temp.folder, temp.folder, ref.fasta)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'SetNmMdAndUqTags'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/fix.bam

   if (keep.unBQSR) {
      from.str = sprintf('%s/fix.bam', temp.folder)
      to.str = sprintf('%s/%s.unBQSR.bam', output.folder, sample.name.str)
      file.copy(from = from.str, to = to.str)

      from.str = sprintf('%s/fix.bai', temp.folder)
      to.str = sprintf('%s/%s.unBQSR.bai', output.folder, sample.name.str)
      file.copy(from = from.str, to = to.str)
   }
   # output: {output.folder}/{sample.name.str}.unBQSR.bam

   # ==========================BaseRecalibrator  output: {temp.folder}/0001.report ... {temp.folder}/0017.report==========================
   # remove all existed report files
   report.files.c = list.files(temp.folder, pattern = '.report$', full.names = T)
   unlink(report.files.c, force = T)  # important

   group.list = CreateSequenceGroupingTSV(ref.dict)
   input.list = list()
   for (i in 1:length(group.list$mapped)) {
      input.list[[i]] = c(no = sprintf('%04d', i), interval = group.list$mapped[i])
   }
   # input.list[[1]]  = c(no = "0001", interval = "-L 1")
   # ......
   # input.list[[17]] = c(no = "0017", interval = "-L 19 -L 22 -L 21 -L GL000192.1 -L GL000225.1 -L GL000194.1 .....")

   # named.vector = c(no = '0012', interval = '-L 18 -L 20 -L Y')
   ..BaseRecalibrator <- function(named.vector) {
      prefix = named.vector['no']  # '0012'

      input.interval.str = named.vector['interval']  # '-L 18 -L 20 -L Y'
      input.sites.str = paste0(' --known-sites ', known.sites.VCFs, collapse = '')
      command.str = sprintf('%s --java-options -Xms8G BaseRecalibrator -R %s -I %s/fix.bam -O %s/%s.report %s %s --use-original-qualities', gatk, ref.fasta, temp.folder, temp.folder, prefix, input.sites.str, input.interval.str)
      message(sprintf('\n[%s] %s: %s - %s', Sys.time(), sample.name.str, 'BaseRecalibrator', prefix))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
      # output:  {temp.folder}/0001.report
   }

   r = parallel::mclapply(input.list, FUN = ..BaseRecalibrator, mc.cores = scatter.count)
   if (all(r == 0)) {
      message(sprintf('\n[%s] SUCCESS: %s BaseRecalibrator', Sys.time(), sample.name.str))
   } else {
      message(sample.name.str, ' BaseRecalibrator:')
      print(r)
      return(1)
   }
   # output: {temp.folder}/0001.report ........ {temp.folder}/0017.report

   # ==========================GatherBqsrReports output : {temp.folder}/BQSR.recal.data.csv==========================
   report.files.c = list.files(temp.folder, pattern = '.report$', full.names = T)
   if (length(report.files.c) == 0) {
      message('No Bqsr Reports Found.')
      quit(save = 'no', status = 1)
   } else {
      message(sprintf('Gather %s reports.', length(report.files.c)))
   }

   input.report.c = paste0(' -I ', report.files.c, collapse = '')
   command.str = sprintf('%s --java-options -Xms8G GatherBQSRReports %s -O %s/BQSR.recal.data.csv', gatk, input.report.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'GatherBQSRReports'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output : {temp.folder}/BQSR.recal.data.csv

   # ==========================ApplyBQSR         output: {temp.folder}/0001.recalibrated.bam  ... {temp.folder}/0017.recalibrated.bam==========================
   # named.vector = c(no = '0012', interval = '-L 18 -L 20 -L Y')
   ..ApplyBQSR <- function(named.vector) {
      prefix = named.vector['no']  # '0012'

      input.interval.str = named.vector['interval']  # '-L 18 -L 20 -L Y'
      command.str = sprintf('%s --java-options -Xms8G ApplyBQSR -R %s -I %s/fix.bam -bqsr %s/BQSR.recal.data.csv -O %s/%s.recalibrated.bam %s --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --use-original-qualities --emit-original-quals true', gatk, ref.fasta, temp.folder, temp.folder, temp.folder, prefix, input.interval.str)
      message(sprintf('\n[%s] %s: %s - %s', Sys.time(), sample.name.str, 'ApplyBQSR', prefix))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }

   }
   # output : {temp.folder}/0012.recalibrated.bam

   # remove all existed recalibrated bam files
   recalibrated.bam.c = list.files(temp.folder, pattern = '.recalibrated.bam$', full.names = T)
   unlink(recalibrated.bam.c, force = T)

   input.list[[length(input.list) + 1]] = c(no = sprintf('%04d', length(input.list) + 1), interval = '-L unmapped')

   r = parallel::mclapply(input.list, FUN = ..ApplyBQSR, mc.cores = scatter.count)
   if (all(r == 0)) {
      message(sprintf('\n[%s] SUCCESS: %s ApplyBQSR', Sys.time(), sample.name.str))
   } else {
      message(sample.name.str, ' ApplyBQSR: ', r)
      return(1)
   }
   # output: {temp.folder}/0001.recalibrated.bam  ... {temp.folder}/0018.recalibrated.bam

   # ==========================GatherBamFiles    output: {output.folder}/{sample.name.str}.bam==========================
   recalibrated.bam.c = list.files(temp.folder, pattern = '.recalibrated.bam$', full.names = T)
   if (length(recalibrated.bam.c) == 0) {
      message('No Recalibrated bam files Found.')
      quit(save = 'no', status = 1)
   }

   input.bam.c = paste0(' -I ', recalibrated.bam.c, collapse = '')
   command.str = sprintf('%s --java-options -Xms8G GatherBamFiles %s -O %s/recalibrated.bam --CREATE_INDEX false', gatk, input.bam.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'GatherBamFiles'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output:  reports : {temp.folder}/recalibrated.bam

   command.str = sprintf('%s --java-options -Xms8G SortSam --SORT_ORDER coordinate -I %s/recalibrated.bam -O %s/%s.bam --CREATE_INDEX true', gatk, temp.folder, output.folder, sample.name.str)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'output'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      unlink(temp.folder, recursive = T)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {output.folder}/{sample.name.str}.bam

   return(0)
}


# run in parallele mode
# temp.folder = sprintf('temp/somatic.snp/%s', sample.name.str)
# named.vector = c(tumor.bam = '~/projects/test/bamfiles/1801167-PT.bam', normal.bam = '~/projects/test/bamfiles/1801167-PN.bam') normal.bam could be NA
# gatk = '~/bin/gatk-4.1.8.1/gatk'
# intervals = c('temp/somatic.snp/1801169-PT/0000-scattered.interval_list' .... '0007-scattered.interval_list')
# scatter.count = 8
# ref.fasta = '~/db/mutect2_support/b37/human_g1k_v37.fasta'
# ref.dict = '~/db/mutect2_support/b37/human_g1k_v37.dict'
# variants.for.contamination = '~/db/mutect2_support/small_exac_common_3_b37.vcf.gz'
# realignment.index.bundle = '~/db/mutect2_support/realignment_index_bundle_hg19.img'
# pon = '~/projects/EHCC/Run_Mutect2_PON/pon.vcf'
# gnomad = '~/db/mutect2_support/af-only-gnomad.raw.sites.b37.vcf.gz'
# gga.vcf = ''
# realignment.extra.args = ''
# m2.extra.args = ''
# getpileupsummaries.extra.args = ''
# filter.mutect2.extra.args = ''
# r = ..Mutect2(named.vector = named.vector, gatk = gatk, interval = interval, scatter.count = scatter.count, ref.fasta = ref.fasta, ref.dict = ref.dict, variants.for.contamination = variants.for.contamination, pon = pon, gnomad = gnomad, realignment.index.bundle = realignment.index.bundle)
# temp.folder: "temp/somatic.snp/1801169-PT"
# output: {output.folder}/{sample.name.str}.vcf.gz
# return 0
..Mutect2 <- function(named.vector, gatk, gatk.safe, intervals, scatter.count, ref.fasta, ref.dict, variants.for.contamination, pon = '', gnomad = '', gga.vcf = '', realignment.index.bundle = '', m2.extra.args = '', getpileupsummaries.extra.args = '', filter.mutect2.extra.args = '', realignment.extra.args = '') {
   Sys.sleep(runif(1, 0.1, 2))

   tumor.bam = normalizePath(named.vector['tumor.bam'])
   normal.bam = ifelse(is.na(named.vector['normal.bam']), yes = NA, no = normalizePath(named.vector['normal.bam'])) # could be NA

   if (is.na(tumor.bam)) {
      message('No tumor for this item. Skip')
      return('No tumor for this item')
   }

   gatk = gatk
   gatk.safe = gatk.safe
   interval.files.c = intervals
   scatter.count = scatter.count
   ref.fasta = ref.fasta
   ref.dict = ref.dict
   variants.for.contamination = variants.for.contamination
   realignment.index.bundle = realignment.index.bundle

   pon = pon
   gnomad = gnomad
   gga.vcf = gga.vcf
   m2.extra.args = m2.extra.args
   getpileupsummaries.extra.args = getpileupsummaries.extra.args
   filter.mutect2.extra.args = filter.mutect2.extra.args

   sample.name.str = tools::file_path_sans_ext(basename(tumor.bam))  # 1801169-PT

   temp.folder = sprintf('temp/somatic.snp/%s', sample.name.str) # "temp/somatic.snp/1801169-PT"
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = './somatic.snp'
   dir.create(output.folder, recursive = T, showWarnings = F)

   # ============M2                        output: {temp.folder}/NNNN/bamout.bam|f1r2.tar.gz|output.vcf|output.vcf.stat|tumor-pileups.tablenormal-pileups.table(optional)==========================
   # interval = 'temp/somatic.snp/NNNN-scattered.interval_list'
   ..M2 <- function(interval) {

      prefix = str_extract(basename(interval), pattern = '\\d+')  # NNNN
      temp.folder = sprintf('%s/%s', temp.folder, prefix)  # "temp/somatic.snp/1801169-PT/NNNN"
      dir.create(temp.folder, recursive = T, showWarnings = F)

      # Mutect2 output: {temp.folder}/NNNN/bamout.bam|f1r2.tar.gz|output.vcf|output.vcf.stats=========================
      tumor.command.line = sprintf('-I %s -tumor %s', tumor.bam, sample.name.str)
      if (!is.na(normal.bam)) {
         normal.command.line = sprintf('-I %s -normal %s', normal.bam, file_path_sans_ext(basename(normal.bam)))
      } else {
         normal.command.line = ''
      }

      command.str = sprintf('%s --java-options -Xmx8g Mutect2 -R %s %s %s -L %s --bam-output %s/bamout.bam --f1r2-tar-gz %s/f1r2.tar.gz %s %s %s -O %s/output.vcf %s',
                            gatk,
                            ref.fasta,
                            tumor.command.line,
                            normal.command.line,
                            interval,
                            temp.folder,
                            temp.folder,
                            ifelse(pon == '', '', sprintf('-pon %s', pon)),
                            ifelse(gnomad == '', '', sprintf('--germline-resource %s', gnomad)),
                            ifelse(gga.vcf == '', '', sprintf('--alleles %s', gga.vcf)),
                            temp.folder,
                            m2.extra.args)
      message(sprintf('\n[%s] %s - %s: %s', Sys.time(), sample.name.str, prefix, 'Mutect2'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }

      # output: 1, temp/1801167-PT/NNNN/bamout.bam
      #         2, temp/1801167-PT/NNNN/f1r2.tar.gz
      #         3, temp/1801167-PT/NNNN/output.vcf
      #         4, temp/1801167-PT/NNNN/output.vcf.stats

      # GetPileupSummaries output: {temp.folder}/NNNN/tumor-pileups.table|normal-pileups.table (optional)=========================

      # for tumor
      command.str = sprintf('%s --java-options -Xmx8g GetPileupSummaries -R %s -I %s --interval-set-rule INTERSECTION -L %s -L %s -V %s -O %s/tumor-pileups.table %s', gatk, ref.fasta, tumor.bam, interval, variants.for.contamination, variants.for.contamination, temp.folder, getpileupsummaries.extra.args)
      message(sprintf('\n[%s] %s - %s: %s', Sys.time(), sample.name.str, prefix, 'GetPileupSummaries(tumor)'))

      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
      # output : temp/1801167-PT/NNNN/tumor-pileups.table

      # for normal
      if (!is.na(normal.bam)) {
         command.str = sprintf('%s --java-options -Xmx8000m GetPileupSummaries -R %s -I %s --interval-set-rule INTERSECTION -L %s -L %s -V %s -O %s/normal-pileups.table %s', gatk, ref.fasta, normal.bam, interval, variants.for.contamination, variants.for.contamination, temp.folder, getpileupsummaries.extra.args)
         message()
         message(sample.name.str, ': GetPileupSummaries(normal) - ', basename(interval))
         r = system(command.str, ignore.stderr = T, ignore.stdout = T)
         if (r == 0) {
            message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         } else {
            message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
            return(1)
         }
         # output : temp/1801167-PT/NNNN/normal-pileups.table
      }

      tumor.table.str = sprintf('%s/tumor-pileups.table', temp.folder)
      if (file.exists(tumor.table.str) & file.size(tumor.table.str) != 0) {
         return(0)
      } else {
         return(1)
      }
   } # M2 ends here

   r = parallel::mclapply(interval.files.c, FUN = ..M2, mc.cores = scatter.count)
   message()
   print(paste0(interval.files.c, '(', r, ')'))
   if (all(r == 0)) {
      message(sprintf('\n[%s] SUCCESS: %s M2', Sys.time(), sample.name.str))
   } else {
      print(r)
      message(sample.name.str, ' M2: ', r)
      return(1)
   }
   # output: 1, {temp.folder}/NNNN/bamout.bam   .....
   #         2, {temp.folder}/NNNN/f1r2.tar.gz  .....
   #         3, {temp.folder}/NNNN/output.vcf   .....
   #         4, {temp.folder}/NNNN/output.vcf.stats
   #         5, {temp.folder}/NNNN/tumor-pileups.table  .....
   #         6, {temp.folder}/NNNN/normal-pileups.table (optional)   .....


   #temp.folder = sprintf('temp/somatic.snp/%s', sample.name.str) # '1801169-PT'
   prefix.c = str_extract(basename(interval.files.c), pattern = '\\d+')  # c('0001', '0002' ......)
   # ============MergeReadOrientationModel output: {temp.folder}/artifact-priors.tar.gz==========================
   files.c = paste0(temp.folder, '/', prefix.c, '/f1r2.tar.gz')
   input.c = str_trim(paste0(' -I ', files.c, collapse = ''))
   command.str = sprintf('%s --java-options -Xmx8g LearnReadOrientationModel %s -O %s/artifact-priors.tar.gz', gatk, input.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'LearnReadOrientationModel (Merge f1r2.tar.gz ...)'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/artifact-priors.tar.gz

   # ============MergeVCFs                 output: {temp.folder}/output.vcf==========================
   files.c = paste0(temp.folder, '/', prefix.c, '/output.vcf')
   input.c = str_trim(paste0(' -I ', files.c, collapse = ''))

   command.str = sprintf('%s --java-options -Xmx8g MergeVcfs %s -O %s/output.vcf', gatk, input.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MergeVcfs'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/output.vcf

   # ============MergeBamOuts              output: {temp.folder}/bam.out.bam==========================
   files.c = paste0(temp.folder, '/', prefix.c, '/bamout.bam')
   input.c = str_trim(paste0(' -I ', files.c, collapse = ''))

   command.str = sprintf('%s --java-options -Xmx8g GatherBamFiles -R %s %s -O %s/unsorted.out.bam', gatk, ref.fasta, input.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MergeBamOuts(merge)'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/unsorted.out.bam

   command.str = sprintf('%s --java-options -Xmx8g SortSam -I %s/unsorted.out.bam -O %s/bam.out.bam --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT --CREATE_INDEX', gatk, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MergeBamOuts(sort)'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/bam.out.bam

   # ============MergeStats                output: {temp.folder}/merged.stats==========================
   files.c = paste0(temp.folder, '/', prefix.c, '/output.vcf.stats')
   input.c = str_trim(paste0(' -stats ', files.c, collapse = ''))

   command.str = sprintf('%s --java-options -Xmx8g MergeMutectStats %s -O %s/merged.stats', gatk, input.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MergeMutectStats'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/merged.stats

   # ============MergePileupSummaries      output: {temp.folder}/MergeTumorPileups.table|MergeNormalPileups.table(optional)==========================
   files.c = paste0(temp.folder, '/', prefix.c, '/tumor-pileups.table')
   input.c = str_trim(paste0(' -I ', files.c, collapse = ''))

   command.str = sprintf('%s --java-options -Xmx4g GatherPileupSummaries --sequence-dictionary %s %s -O %s/MergeTumorPileups.table', gatk, ref.dict, input.c, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MergePileupSummaries(tumor)'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output : {temp.folder}/MergeTumorPileups.table

   if (!is.na(normal.bam)) {
      files.c = paste0(temp.folder, '/', prefix.c, '/normal-pileups.table')
      input.c = str_trim(paste0(' -I ', files.c, collapse = ''))

      command.str = sprintf('%s --java-options -Xmx4g GatherPileupSummaries --sequence-dictionary %s %s -O %s/MergeNormalPileups.table', gatk, ref.dict, input.c, temp.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'MergePileupSummaries(normal)'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
      # output : {temp.folder}/MergeNormalPileups.table
   }

   # ============CalculateContamination    output: {temp.folder}/contamination.table|segments.table==========================
   normal.merged.table.str = sprintf('%s/MergeNormalPileups.table', temp.folder)

   input.match.str = ifelse(!is.na(named.vector['normal.bam']) & file.exists(normal.merged.table.str), yes = sprintf('-matched %s', normal.merged.table.str), no = '')
   command.str = sprintf('%s --java-options -Xmx8g CalculateContamination -I %s/MergeTumorPileups.table -O %s/contamination.table --tumor-segmentation %s/segments.table %s', gatk, temp.folder, temp.folder, temp.folder, input.match.str)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CalculateContamination'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: 1, {temp.folder}/contamination.table
   #         2, {temp.folder}/segments.table

   # ============Filter                    output: {temp.folder}/filtered.vcf|filtering.stats==========================
   command.str = sprintf('%s --java-options -Xmx8g FilterMutectCalls -R %s -V %s/output.vcf -O %s/filtered.vcf --contamination-table %s/contamination.table --tumor-segmentation %s/segments.table --ob-priors %s/artifact-priors.tar.gz -stats %s/merged.stats --filtering-stats %s/filtering.stats %s', gatk, ref.fasta, temp.folder, temp.folder, temp.folder, temp.folder, temp.folder, temp.folder, temp.folder, filter.mutect2.extra.args)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'FilterMutectCalls'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output:1, {temp.folder}/filtered.vcf
   #        2, {temp.folder}/filtering.stats

   # ============FilterAlignmentArtifacts  output: {temp.folder}/final.vcf.gz==========================
   command.str = sprintf('%s --java-options -Xmx32g FilterAlignmentArtifacts -R %s -V %s/filtered.vcf -I %s/bam.out.bam -O %s/final.vcf --bwa-mem-index-image %s --create-output-variant-index false %s', gatk, ref.fasta, temp.folder, temp.folder, temp.folder, realignment.index.bundle, realignment.extra.args)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'FilterAlignmentArtifacts'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      message(sprintf('\nTry GATK fail-safe version (%s).', gatk.safe))
      command.str = sprintf('%s --java-options -Xmx32g FilterAlignmentArtifacts -R %s -V %s/filtered.vcf -I %s/bam.out.bam -O %s/final.vcf --bwa-mem-index-image %s --create-output-variant-index false %s', gatk.safe, ref.fasta, temp.folder, temp.folder, temp.folder, realignment.index.bundle, realignment.extra.args)
      message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'FilterAlignmentArtifacts (Fail safe)'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         message('\nGATK fail-safe version does not work.')
         return(1)
      }
   }
   # output: {temp.folder}/final.vcf.gz

   # ============SortVcf                   output: {output.folder}/{sample.name.str}.vcf.gz==========================
   command.str = sprintf('%s SortVcf -I %s/final.vcf -O %s/%s.vcf.gz', gatk, temp.folder, output.folder, sample.name.str)
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'Sort the final output'))
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      return(0)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   #output: {output.folder}/{sample.name.str}.vcf.gz
}


# run in parallele mode
# bam.file = '~/projects/test/bamfiles/1801167-PN.bam'
# gatk = '~/bin/gatk-4.1.8.1/gatk'
# interval = '~/db/BEDs/Agilent_SureSelect_Human_All_Exon_V6_hg19_noXY_noChr.interval_list'
# ref.fasta = '~/db/mutect2_support/b37/human_g1k_v37.fasta'
# output  temp/variant_germline/1801169-PN/HaplotypeCaller.g.vcf.gz
# temp/variant_germline/1801167-PN
# return 0
..HaplotypeCallerGvcf <- function(bam.file, gatk, interval, ref.fasta, str.table) {
   Sys.sleep(runif(1, 0.1, 2))

   sample.name.str = tools::file_path_sans_ext(basename(bam.file)) # "1801167-PN"

   temp.folder = sprintf('temp/germline.snp/%s', sample.name.str)
   dir.create(temp.folder, recursive = T, showWarnings = F)


   # =========================CalibrateDragstrModel    output: {temp.folder}/STR.txt=========================
   # ~/bin/gatk-4.2.0.0/gatk CalibrateDragstrModel -I bamfiles/962290-DT.bam -O str.test.txt -R ~/db/mutect2_support/b37/human_g1k_v37.fasta -str ~/db/mutect2_support/b37/human_g1k_v37.STR.table.zip
   command.str = sprintf('%s CalibrateDragstrModel -I %s -O %s/STR.txt -R %s -str %s', gatk, bam.file, temp.folder, ref.fasta, str.table)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CalibrateDragstrModel'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/STR.txt


   # =========================HaplotypeCaller          output: {temp.folder}/HaplotypeCaller.g.vcf.gz=========================
   contamination.file.str = sprintf('temp/somatic.snp/%s/contamination.table', sample.name.str)
   if (file.exists(contamination.file.str)) {
      contamination.float = read.csv(contamination.file.str, sep = '\t')$contamination
   } else {
      contamination.float = 0
   }

   command.str = sprintf('%s --java-options -Xmx4g HaplotypeCaller -R %s -I %s -L %s -O %s/HaplotypeCaller.g.vcf.gz -contamination %s --dragstr-params-path %s/STR.txt -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF --dragen-mode true', gatk, ref.fasta, bam.file, interval, temp.folder, contamination.float, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'HaplotypeCaller'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      return(0)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/HaplotypeCaller.g.vcf.gz
}
# output: temp/variant_germline/1801169-PN/HaplotypeCaller.g.vcf.gz


# Running this function in parallel mode is NOT allowed
# gVCFs.files.c = c('temp/1801169-PT/HaplotypeCaller.g.vcf.gz', 'temp/1801169-PN/HaplotypeCaller.g.vcf.gz', 'temp/1801167-PT/HaplotypeCaller.g.vcf.gz', 'temp/1801167-PN/HaplotypeCaller.g.vcf.gz')
# gatk = '~/bin/gatk-4.1.8.1/gatk'
# interval = '~/db/BEDs/Agilent_SureSelect_Human_All_Exon_V6_hg19_noXY_noChr.interval_list'
# scatter.count = 16
# ref.fasta = '~/db/mutect2_support/b37/human_g1k_v37.fasta'
# ref.dict = '~/db/mutect2_support/b37/human_g1k_v37.dict'
# dbsnp = 'dbsnp_138.b37.vcf.gz'
# mills.resource = 'Mills_and_1000G_gold_standard.indels.b37.sites.vcf'
# axiomPoly.resource = 'Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz'
# hapmap.resource = 'hapmap_3.3.b37.vcf.gz'
# omni.resource = '1000G_omni2.5.b37.vcf.gz'
# one.thousand.genomes.resource ='1000G_phase1.snps.high_confidence.b37.vcf.gz'
# suitable when samples are no greater than 1000
# temp.folder: temp/germline.snp
# output ./germline.snp/jointGenotyping.vcf
# return 0
..JointGenotyping <- function(gVCFs.files.c, gatk, gatk.safe, interval, scatter.count, sample.each.time, ref.fasta, dbsnp, mills.resource, axiomPoly.resource, hapmap.resource, omni.resource, one.thousand.genomes.resource, max.gaussians.indel = 4, max.gaussians.snp = 6) {

   gVCFs.files.c = gVCFs.files.c
   gatk = gatk
   gatk.safe = gatk.safe
   interval = interval
   scatter.count = scatter.count
   ref.fasta = ref.fasta

   dbsnp = dbsnp
   mills.resource = mills.resource
   axiomPoly.resource = axiomPoly.resource
   hapmap.resource = hapmap.resource
   omni.resource = omni.resource
   one.thousand.genomes.resource = one.thousand.genomes.resource

   temp.folder = 'temp/germline.snp'
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = 'germline.snp'
   dir.create(output.folder, recursive = T, showWarnings = F)

   #==============SplitIntervals          output: {temp.folder}/NNNN-scattered.interval_list===============
   unpadded_intervals.c = list.files(temp.folder, '^\\d+-scattered\\.interval_list$', full.names = T)
   unlink(unpadded_intervals.c)

   command.str = sprintf('%s --java-options -Xms4g SplitIntervals -scatter %s -L %s -O %s/ -R %s -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --interval-merging-rule OVERLAPPING_ONLY', gatk, scatter.count, interval, temp.folder, ref.fasta)
   message(sprintf('\n[%s] %s', Sys.time(), 'SplitIntervals'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit(save = 'no', status = r)
   }
   # output: temp/variant_germline/0000-scattered.interval_list, ..... temp/variant_germline/0008-scattered.interval_list
   unpadded_intervals.c = list.files(temp.folder, '^\\d+-scattered\\.interval_list$', full.names = T)
   interval.list = list()
   for (file in unpadded_intervals.c) {
      interval.list[[basename(file)]] = file
   }

   #==============GenotypeAndHardFilter   output: {temp.folder}/NNNN/filtered.vcf|sites.only.vcf.gz===========
   # subinterval = 'temp/variant_germline/0000-scattered.interval_list'
   ..GenotypeAndHardFilter <- function(subinterval) {
      Sys.sleep(runif(1, 0.1, 2))

      prefix = str_extract(basename(subinterval), '\\d+')  # 'NNNN'
      message(sprintf('\nDeleting %s/genomicsdb_%s ...', temp.folder, prefix))
      unlink(sprintf('%s/genomicsdb_%s', temp.folder, prefix), recursive = T)  # important

      temp.interval.folder = sprintf('%s/%s', temp.folder, prefix)  # temp/germline.snp/NNNN
      dir.create(temp.interval.folder, recursive = T, showWarnings = F)

      # GenomicsDBImport
      dir.create(sprintf('%s/GenomicsDBImport.temp/%s', temp.folder, prefix), recursive = T, showWarnings = F)

      input.c = paste0(' -V ', gVCFs.files.c, collapse = '')
      command.str = sprintf('%s --java-options -Xms8g GenomicsDBImport --tmp-dir %s/GenomicsDBImport.temp/%s --genomicsdb-workspace-path %s/genomicsdb_%s --batch-size 50 -L %s %s --reader-threads %s --merge-input-intervals --consolidate', gatk, temp.folder, prefix, temp.folder, prefix, subinterval, input.c, sample.each.time)
      message(sprintf('\n[%s] %s: %s', Sys.time(), 'JointGenotyping GenomicsDBImport', prefix))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         unlink(sprintf('%s/GenomicsDBImport.temp/%s', temp.folder, prefix))
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }
      # output: {temp.folder}/genomicsdb_0000

      # GenotypeGVCFs
      command.str = sprintf('%s --java-options -Xms8g GenotypeGVCFs -R %s -O %s/raw.vcf -V gendb://%s/genomicsdb_%s -L %s -D %s -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals --merge-input-intervals', gatk, ref.fasta, temp.interval.folder, temp.folder, prefix, subinterval, dbsnp)
      message(sprintf('\n[%s] %s: %s', Sys.time(), 'JointGenotyping GenotypeGVCFs', prefix))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         unlink(sprintf('%s/genomicsdb_%s', temp.folder, prefix), recursive = T)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }
      # output: {temp.interval.folder}/raw.vcf

      # VariantFiltration
      command.str = sprintf('%s --java-options -Xms4g VariantFiltration -V %s/raw.vcf -O %s/filtered.vcf --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet', gatk, temp.interval.folder, temp.interval.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), 'JointGenotyping VariantFiltration', prefix))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }
      # output: {temp.interval.folder}/filtered.vcf

      # MakeSitesOnlyVcf
      command.str = sprintf('%s --java-options -Xms4g MakeSitesOnlyVcf -I %s/filtered.vcf -O %s/sites.only.vcf.gz', gatk, temp.interval.folder, temp.interval.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), 'JointGenotyping MakeSitesOnlyVcf', prefix))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }
      # output: {temp.interval.folder}/sites.only.vcf.gz

      return(0)
   }
   # output: 1, temp/variant_germline/0000/filtered.vcf
   #         2, temp/variant_germline/0000/sites.only.vcf.gz

   r = parallel::mclapply(interval.list, FUN = ..GenotypeAndHardFilter, mc.cores = scatter.count)
   message()
   print(r)
   message(r)
   # output: temp/variant_germline/0000/sites.only.vcf   temp/variant_germline/000X/sites.only.vcf
   # output: temp/variant_germline/0000/filtered.vcf   temp/variant_germline/000X/filtered.vcf

   #==============SitesOnlyGatherVcf      output: {temp.folder}/sites.only.vcf.gz===========
   prefix.c = str_extract(basename(unpadded_intervals.c), '\\d+')   # '0000', '0001', ... '000N'
   sites.only.vcf.c = paste0(temp.folder, '/', prefix.c, '/sites.only.vcf.gz')
   input.c = paste0(' -I ', sites.only.vcf.c, collapse = '')

   command.str = sprintf('%s --java-options -Xms6g GatherVcfsCloud %s --ignore-safety-checks --gather-type BLOCK -O %s/temp.vcf.gz', gatk, input.c, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'GatherVcfsCloud'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit('no', status = r)
   }

   command.str = sprintf('%s --java-options -Xms6g SortVcf -I %s/temp.vcf.gz -O %s/sites.only.vcf.gz', gatk, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'SortVcf'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit('no', status = r)
   }
   # output: {temp.folder}/sites.only.vcf.gz

   #==============IndelsVariantRecalibrator & SNPsVariantRecalibrator   output: {temp.folder}/recaliberation.indel.vcf|tranches.indel.txt|recaliberation.snp.vcf|tranches.snp.txt===========

   command.list = list()

   # -mode INDEL
   command.str = sprintf('%s --java-options -Xms24g VariantRecalibrator -V %s/sites.only.vcf.gz -O %s/recaliberation.indel.vcf --tranches-file %s/tranches.indel.txt --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an DP -an FS -an MQRankSum -an QD -an ReadPosRankSum -an SOR -mode INDEL --max-gaussians %s -resource:mills,known=false,training=true,truth=true,prior=12 %s -resource:axiomPoly,known=false,training=true,truth=false,prior=10 %s -resource:dbsnp,known=true,training=false,truth=false,prior=2 %s --use-allele-specific-annotations', gatk, temp.folder, temp.folder, temp.folder, max.gaussians.indel, mills.resource, axiomPoly.resource, dbsnp)

   command.2.str = sprintf('%s --java-options -Xms24g VariantRecalibrator -V %s/sites.only.vcf.gz -O %s/recaliberation.indel.vcf --tranches-file %s/tranches.indel.txt --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an DP -an FS -an MQRankSum -an QD -an ReadPosRankSum -an SOR -mode INDEL --max-gaussians %s -resource:mills,known=false,training=true,truth=true,prior=12 %s -resource:axiomPoly,known=false,training=true,truth=false,prior=10 %s -resource:dbsnp,known=true,training=false,truth=false,prior=2 %s --use-allele-specific-annotations', gatk.safe, temp.folder, temp.folder, temp.folder, max.gaussians.indel, mills.resource, axiomPoly.resource, dbsnp)

   command.list$VariantRecalibrator.INDEL = c(prefer = command.str, fail.safe = command.2.str)

   # -mode SNP
   command.str = sprintf('%s --java-options -Xms24g VariantRecalibrator -V %s/sites.only.vcf.gz -O %s/recaliberation.snp.vcf --tranches-file %s/tranches.snp.txt --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an DP -an FS -an MQ -an MQRankSum -an QD -an ReadPosRankSum -an SOR -mode SNP --max-gaussians %s -resource:hapmap,known=false,training=true,truth=true,prior=15 %s -resource:omni,known=false,training=true,truth=true,prior=12 %s -resource:1000G,known=false,training=true,truth=false,prior=10 %s -resource:dbsnp,known=true,training=false,truth=false,prior=7 %s --use-allele-specific-annotations', gatk, temp.folder, temp.folder, temp.folder, max.gaussians.snp, hapmap.resource, omni.resource, one.thousand.genomes.resource, dbsnp)

   command.2.str = sprintf('%s --java-options -Xms24g VariantRecalibrator -V %s/sites.only.vcf.gz -O %s/recaliberation.snp.vcf --tranches-file %s/tranches.snp.txt --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an DP -an FS -an MQ -an MQRankSum -an QD -an ReadPosRankSum -an SOR -mode SNP --max-gaussians %s -resource:hapmap,known=false,training=true,truth=true,prior=15 %s -resource:omni,known=false,training=true,truth=true,prior=12 %s -resource:1000G,known=false,training=true,truth=false,prior=10 %s -resource:dbsnp,known=true,training=false,truth=false,prior=7 %s --use-allele-specific-annotations', gatk.safe, temp.folder, temp.folder, temp.folder, max.gaussians.snp, hapmap.resource, omni.resource, one.thousand.genomes.resource, dbsnp)

   command.list$VariantRecalibrator.SNP = c(prefer = command.str, fail.safe = command.2.str)

   message = sprintf('\n[%s] %s', Sys.time(), 'Germline Variant Recalibrator')

   r = parallel::mclapply(command.list, FUN = run, message = message, quit.on.fail = T, mc.cores = 2)
   message()
   print(r)
   if (all(r == 0)) {
      message(sprintf('\n[%s] SUCCESS: Germline Variant Recalibrator', Sys.time()))
   } else {
      print(r)
   }
   # output: 1, {temp.folder}/recaliberation.indel.vcf
   #         2, {temp.folder}/tranches.indel.txt
   #         3, {temp.folder}/recaliberation.snp.vcf
   #         4, {temp.folder}/tranches.snp.txt


   #==============ApplyRecalibration      output: {temp.folder}/NNNN/snp.recalibrated.vcf===========
   prefix.c = str_extract(basename(unpadded_intervals.c), '\\d+')   # '0000', '0001', ... 'NNNN'
   vcf.files.c = normalizePath(paste0(temp.folder, '/', prefix.c, '/filtered.vcf'))
   # "{temp.folder}/0000/filtered.vcf"
   # "{temp.folder}/0001/filtered.vcf"
   #  ....
   # "{temp.folder}/NNNN/filtered.vcf"
   print(vcf.files.c)
   message(sprintf('\n Apply recalibration to %s vcf files', length(vcf.files.c)))

   ..ApplyRecalibration <- function(vcf) {
      Sys.sleep(runif(1, 0.1, 2))

      # ApplyVQSR INDEL
      command.str = sprintf('%s --java-options -Xms5g ApplyVQSR -O %s/indel.recalibrated.vcf -V %s --recal-file %s/recaliberation.indel.vcf --tranches-file %s/tranches.indel.txt --use-allele-specific-annotations --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL', gatk, dirname(vcf), vcf, temp.folder, temp.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), 'ApplyVQSR INDEL', vcf))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
      # output : {temp.folder}/NNNN/indel.recalibrated.vcf


      # ApplyVQSR SNP
      command.str = sprintf('%s --java-options -Xms5g ApplyVQSR -O %s/snp.recalibrated.vcf -V %s/indel.recalibrated.vcf --recal-file %s/recaliberation.snp.vcf --tranches-file %s/tranches.snp.txt --use-allele-specific-annotations --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode SNP', gatk, dirname(vcf), dirname(vcf), temp.folder, temp.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), 'ApplyVQSR SNP', vcf))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
      # output: {temp.folder}/NNNN/snp.recalibrated.vcf
   }

   r = parallel::mclapply(vcf.files.c, FUN = ..ApplyRecalibration)
   message('\nApplyRecalibration reuslt: \n')
   print(paste0(vcf.files.c, ' (', r, ')'))
   message(r)
   # output: {temp.folder}/0000/snp.recalibrated.vcf ... {temp.folder}/NNNN/snp.recalibrated.vcf

   #==============GatherVcfsCloud         output: {output.folder}/jointGenotyping.vcf===========
   vcf.files.c = normalizePath(paste0(temp.folder, '/', prefix.c, '/snp.recalibrated.vcf'))
   input.c = paste0(' -I ', vcf.files.c, collapse = '')

   command.str = sprintf('%s --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type AUTOMATIC %s -O %s/recalibrated.vcf', gatk, input.c, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'GatherVcfsCloud'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit(save = 'no', status = r)
   }
   # output: {temp.folder}/recalibrated.vcf

   command.str = sprintf('%s SortVcf -I %s/recalibrated.vcf -O %s/jointGenotyping.vcf', gatk, temp.folder, output.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'SortVcf Final'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      return(0)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit(save = 'no', status = r)
   }

}


# Running this function in parallel mode is not allowed
# sample.normal.list =
# gatk = '~/bin/gatk-4.1.8.1/gatk'
# sample.each.time = 8
# ref.fasta = '~/db/mutect2_support/b37/human_g1k_v37.fasta'
# ref.dict = '~/db/mutect2_support/b37/human_g1k_v37.dict'
# interval = '~/db/BEDs/Agilent_SureSelect_Human_All_Exon_V6_hg19_noXY_noChr.interval_list'
# common.sites = '~/db/mutect2_support/small_exac_common_3_b37.vcf.gz'
# blacklist_interval = '~/db/gatk_cnv_support/b37/CNV_and_centromere_blacklist.hg19.list'
# segmental_duplication_track_bed = '~/db/gatk_cnv_support/b37/germline-copy-number_hg19.nochr.SegDups.elements.no_gl.bed'
# mappability_track_bed = '~/db/gatk_cnv_support/b37/germline-copy-number_hg19.nochr.k100.umap.single.merged.bed'
# output: {temp.folder}/cnv.somatic.pon.hdf5
# return 0
..CNV_Somatic_Panel <- function(sample.normal.list, gatk, sample.each.time, ref.fasta, ref.dict, interval, common.sites, blacklist_interval = '', segmental_duplication_track_bed = '', mappability_track_bed = '') {

   gatk = gatk
   sample.each.time = sample.each.time
   ref.fasta = ref.fasta
   ref.dict = ref.dict
   interval = interval
   common.sites = common.sites

   blacklist_interval = blacklist_interval
   segmental_duplication_track_bed = segmental_duplication_track_bed
   mappability_track_bed = mappability_track_bed

   temp.folder = 'temp/somatic.cnv'
   dir.create(temp.folder, recursive = T, showWarnings = F)

   #===================PreprocessIntervals    output: {temp.folder}/preprocessed.interval_list===============
   input.black.str = ifelse(blacklist_interval == '', yes = '', no = sprintf('-XL %s', blacklist_interval))
   command.str = sprintf('%s --java-options -Xmx4g PreprocessIntervals -L %s %s --reference %s --sequence-dictionary %s --output %s/preprocessed.interval_list --padding 250 --bin-length 1000 --interval-merging-rule OVERLAPPING_ONLY', gatk, interval, input.black.str, ref.fasta, ref.dict, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'CNV PON: PreprocessIntervals'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit('no', status = r)
   }
   # output: {temp.folder}/preprocessed.interval_list

   #===================AnnotateIntervals      output: {temp.folder}/annotated.interval_list===============
   input.map.str = ifelse(mappability_track_bed == '', yes = '', no = sprintf('--mappability-track %s', mappability_track_bed))
   input.sd.str = ifelse(segmental_duplication_track_bed == '', yes = '', no = sprintf('--segmental-duplication-track %s', segmental_duplication_track_bed))
   command.str = sprintf('%s --java-options -Xmx8g AnnotateIntervals -R %s -L %s/preprocessed.interval_list -O %s/annotated.interval_list %s %s --feature-query-lookahead 1000000 --interval-merging-rule OVERLAPPING_ONLY', gatk, ref.fasta, temp.folder, temp.folder, input.map.str, input.sd.str)
   message(sprintf('\n[%s] %s', Sys.time(), 'CNV PON: AnnotateIntervals'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit('no', status = r)
   }
   # output: {temp.folder}/annotated.interval_list

   #==================CollectCounts           output: {temp.folder}/{sample.name.str}/counts.hdf5.......===============
   ..CollectCounts <- function(bam.file) {
      Sys.sleep(runif(1, 0.2, 2))

      sample.name.str = file_path_sans_ext(basename(bam.file))
      temp.sample.folder = sprintf('%s/%s', temp.folder, sample.name.str)
      dir.create(temp.sample.folder, recursive = T, showWarnings = F)

      command.str = sprintf('%s --java-options -Xmx4g CollectReadCounts -L %s/preprocessed.interval_list --input %s --reference %s --output %s/counts.hdf5 --format HDF5 --interval-merging-rule OVERLAPPING_ONLY', gatk, temp.folder, bam.file, ref.fasta, temp.sample.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CollectReadCounts PON'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
      # output: {temp.folder}/1801167-PN/counts.hdf5

   }

   r = parallel::mclapply(sample.normal.list, FUN = ..CollectCounts, mc.cores = sample.each.time)
   message()
   print(r)
   message('..CollectCounts: ', r)
   # output: temp/cnv_somatic/1801167-PN/counts.hdf5   ...  temp/cnv_somatic/XXXXXXX-PN/counts.hdf5


   #==================FilterInterval..........output: {temp.folder}/filtered.interval_list         (this step is not for somatic PON creation, but for germline CNV calling)==================
   temp.c = sprintf('--input %s/%s/counts.hdf5', temp.folder, file_path_sans_ext(basename(unlist(sample.normal.list))))
   input.counts.str = paste0(temp.c, collapse = ' ')

   command.str = sprintf('%s --java-options -Xmx8g FilterIntervals -L %s/preprocessed.interval_list -XL %s --annotated-intervals %s/annotated.interval_list --output %s/filtered.interval_list %s --minimum-gc-content 0.1 --maximum-gc-content 0.9 --minimum-mappability 0.9 --maximum-mappability 1.0 --minimum-segmental-duplication-content 0.0 --maximum-segmental-duplication-content 0.5 --low-count-filter-count-threshold 5 --low-count-filter-percentage-of-samples 90.0 --extreme-count-filter-minimum-percentile 1.0 --extreme-count-filter-maximum-percentile 99.0 --extreme-count-filter-percentage-of-samples 90.0 --interval-merging-rule OVERLAPPING_ONLY', gatk, temp.folder, blacklist_interval, temp.folder, temp.folder, input.counts.str)
   message(sprintf('\n[%s] %s', Sys.time(), 'CNV PON: FilterIntervals'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      message('This step is not for somatic PON creation, but for Germline CNV calling. Hence you can still do the somatic CNV calling.')
      message('However this failure will lead to a incapability of Germline CNV calling')
      #quit('no', status = r)
   }
   # output: {temp.folder}/filtered.interval_list


   #==================CollectAllelicCounts    output: {temp.folder}/{sample.name.str}/allelicCounts.tsv (this step is not for somatic PON creation, but for paired somatic CNV calling)===============
   ..CollectAllelicCounts <- function(bam.file) {
      Sys.sleep(runif(1, 0.1, 2))
      sample.name.str = file_path_sans_ext(basename(bam.file))
      temp.sample.folder = sprintf('%s/%s', temp.folder, sample.name.str)
      dir.create(temp.sample.folder, recursive = T, showWarnings = F)

      command.str = sprintf('%s --java-options -Xmx4g CollectAllelicCounts -L %s --input %s --reference %s --output %s/allelicCounts.tsv --minimum-base-quality 20', gatk, common.sites, bam.file, ref.fasta, temp.sample.folder)
      message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CNV PON: CollectAllelicCounts'))

      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
   }

   r = parallel::mclapply(sample.normal.list, FUN = ..CollectAllelicCounts, mc.cores = sample.each.time)
   message()
   print(r)
   message('..CollectAllelicCounts: ', r)
   # output: temp/cnv_somatic/1801167-PN/allelicCounts.tsv  .... temp/cnv_somatic/XXXXXXX-PN/allelicCounts.tsv

   #==================CreateReadCountPanelOfNormals   output: {temp.folder}/cnv.somatic.pon.hdf5===============
   input.count.str = paste0(sprintf(' -I %s/', temp.folder), file_path_sans_ext(basename(unlist(sample.normal.list))), '/counts.hdf5', collapse = '')

   command.str = sprintf('%s --java-options -Xmx8g CreateReadCountPanelOfNormals %s -O %s/cnv.somatic.pon.hdf5 --annotated-intervals %s/annotated.interval_list --minimum-interval-median-percentile 10.0 --maximum-zeros-in-sample-percentage 5.0 --maximum-zeros-in-interval-percentage 5.0 --extreme-sample-median-percentile 2.5 --do-impute-zeros true --extreme-outlier-truncation-percentile 0.1 --number-of-eigensamples 20 --maximum-chunk-size 16777216', gatk, input.count.str, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'CNV PON: CreateReadCountPanelOfNormals'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit('no', status = r)
   }
   # output: {temp.folder}/cnv.somatic.pon.hdf5

   if (file.exists(sprintf('%s/cnv.somatic.pon.hdf5', temp.folder))) {
      return(0)
   } else {
      return(1)
   }

}
# output: temp/somatic.cnv/cnv.somatic.pon.hdf5


# run in parallel mode
# named.vector = c(tumor.bam = '~/projects/test/bamfiles/1801167-DT.bam', normal.bam = '~/projects/test/bamfiles/1801167-DN.bam')
# gatk = '~/bin/gatk-4.1.8.1/gatk'
# sample.each.time = 8
# ref.fasta = '~/db/mutect2_support/b37/human_g1k_v37.fasta'
# ref.dict = '~/db/mutect2_support/b37/human_g1k_v37.dict'
# preprocessed.interval = 'temp/somatic.cnv/preprocessed.interval_list'
# common.sites = '~/db/mutect2_support/small_exac_common_3_b37.vcf.gz'
# return 0
# output: temp/somatic.cnv/{sample.name.str}/modelSegments.called.seg
..CNV_Somatic_Pair_Call <- function(named.vector, gatk, ref.fasta, common.sites, preprocessed.interval) {

   sample.name.str = tools::file_path_sans_ext(basename(named.vector['tumor.bam']))  # 1801167-DT

   if (is.na(named.vector['tumor.bam'])) {
      print(named.vector)
      message('No tumor for this item. Skip')
      return('No tumor for this item. Skip')
   }

   if (is.na(named.vector['normal.bam'])) {
      print(named.vector)
      message('Tumor-only mode (Not recommended).')
   }

   gatk = gatk
   preprocessed.interval = preprocessed.interval
   ref.fasta = ref.fasta
   common.sites = common.sites

   temp.folder = sprintf('temp/somatic.cnv/%s', sample.name.str)  # temp/somatic.cnv/1801167-DT
   dir.create(temp.folder, recursive = T, showWarnings = F)

   #output.folder = sprintf('temp/somatic.cnv/%s', sample.name.str)
   #dir.create(temp.folder, recursive = T, showWarnings = F)

   panel.of.normals = 'temp/somatic.cnv/cnv.somatic.pon.hdf5'

   #===========================CollectReadCounts        output: {temp.folder}/counts.hdf5===========================
   command.str = sprintf('%s --java-options -Xmx4g CollectReadCounts -L %s -I %s -R %s -O %s/counts.hdf5 --format HDF5 --interval-merging-rule OVERLAPPING_ONLY', gatk, preprocessed.interval, named.vector['tumor.bam'], ref.fasta, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CNV somatic: CollectReadCounts'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/counts.hdf5

   #===========================DenoiseReadCounts        output: {temp.folder}/standardizedCR.tsv|denoisedCR.tsv==========================
   command.str = sprintf('%s --java-options -Xmx4g DenoiseReadCounts -I %s/counts.hdf5 --standardized-copy-ratios %s/standardizedCR.tsv --denoised-copy-ratios %s/denoisedCR.tsv --count-panel-of-normals %s', gatk, temp.folder, temp.folder, temp.folder, panel.of.normals)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CNV somatic: DenoiseReadCounts'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }

   # output: 1, temp/cnv_somatic/1801167-PT/standardizedCR.tsv
   #         2, temp/cnv_somatic/1801167-PT/denoisedCR.tsv

   #===========================CollectAllelicCounts     output: {temp.folder}/allelicCounts.tsv===========================
   command.str = sprintf('%s --java-options -Xmx4g CollectAllelicCounts -L %s -I %s -R %s -O %s/allelicCounts.tsv --minimum-base-quality 20', gatk, common.sites, named.vector['tumor.bam'], ref.fasta, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CNV somatic: CollectAllelicCounts'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: 1, {temp.folder}/allelicCounts.tsv



   #===========================ModelSegments            output: {temp.folder}/modelSegments.cr.seg (11 files)===========================
   # if matched normal pair (e.g. temp/cnv_somatic/1448279-DN/allelicCounts.tsv) exists
   if (!is.na(named.vector['normal.bam']) & file.exists(sprintf('temp/cnv_somatic/%s/allelicCounts.tsv', file_path_sans_ext(basename(named.vector['normal.bam']))))) {
      min.total.allele.count.case = 0
      input.normal = sprintf('--normal-allelic-counts temp/cnv_somatic/%s/allelicCounts.tsv --minimum-total-allele-count-normal 30', file_path_sans_ext(basename(named.vector['normal.bam'])))
   } else {
      min.total.allele.count.case = 30
      input.normal = ''
   }

   command.str = sprintf('%s --java-options -Xmx4g ModelSegments --denoised-copy-ratios %s/denoisedCR.tsv --allelic-counts %s/allelicCounts.tsv --output %s %s --output-prefix modelSegments --minimum-total-allele-count-case %s --genotyping-homozygous-log-ratio-threshold -10.0 --genotyping-base-error-rate 0.05 --maximum-number-of-segments-per-chromosome 1000 --kernel-variance-copy-ratio 0.0 --kernel-variance-allele-fraction 0.025 --kernel-scaling-allele-fraction 1.0 --kernel-approximation-dimension 100 --window-size 8 --window-size 16 --window-size 32 --window-size 64 --window-size 128 --window-size 256 --number-of-changepoints-penalty-factor 1.0 --minor-allele-fraction-prior-alpha 25.0 --number-of-samples-copy-ratio 100 --number-of-burn-in-samples-copy-ratio 50 --number-of-samples-allele-fraction 100 --number-of-burn-in-samples-allele-fraction 50 --smoothing-credible-interval-threshold-copy-ratio 2.0 --smoothing-credible-interval-threshold-allele-fraction 2.0 --maximum-number-of-smoothing-iterations 10 --number-of-smoothing-iterations-per-fit 0', gatk, temp.folder, temp.folder, temp.folder, input.normal, min.total.allele.count.case)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CNV somatic: ModelSegments'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   #output: File copy_ratio_only_segments         = "${temp.folder}/modelSegments.cr.seg"
   #        File het_allelic_counts               = "${temp.folder}/modelSegments.hets.tsv"
   #        File normal_het_allelic_counts        = "${temp.folder}/modelSegments.hets.normal.tsv"
   #        File copy_ratio_legacy_segments       = "${temp.folder}/modelSegments.cr.igv.seg"
   #        File allele_fraction_legacy_segments  = "${temp.folder}/modelSegments.af.igv.seg"
   #        File modeled_segments_begin           = "${temp.folder}/modelSegments.modelBegin.seg"
   #        File copy_ratio_parameters_begin      = "${temp.folder}/modelSegments.modelBegin.cr.param"
   #        File allele_fraction_parameters_begin = "${temp.folder}/modelSegments.modelBegin.af.param"
   #        File modeled_segments                 = "${temp.folder}/modelSegments.modelFinal.seg"
   #        File copy_ratio_parameters            = "${temp.folder}/modelSegments.modelFinal.cr.param"
   #        File allele_fraction_parameters       = "${temp.folder}/modelSegments.modelFinal.af.param"


   #===========================CallCopyRatioSegments    output: {temp.folder}/modelSegments.called.seg===========================
   command.str = sprintf('%s --java-options -Xmx4g CallCopyRatioSegments -I %s/modelSegments.cr.seg -O %s/modelSegments.called.seg --neutral-segment-copy-ratio-lower-bound 0.9 --neutral-segment-copy-ratio-upper-bound 1.1 --outlier-neutral-segment-copy-ratio-z-score-threshold 2.0 --calling-copy-ratio-z-score-threshold 2.0', gatk, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'CNV somatic: CallCopyRatioSegments'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # File called_copy_ratio_segments = {temp.folder}/modelSegments.called.seg
   # File called_copy_ratio_legacy_segments = {temp.folder}/modelSegments.called.igv.seg

   CallCopyRatioSegments.str = sprintf('%s/modelSegments.called.seg', temp.folder)
   if (file.exists(CallCopyRatioSegments.str) & file.size(CallCopyRatioSegments.str) != 0) {
      return(0)
   } else {
      return(1)
   }

}
# output: temp/somatic.cnv/{sample.name.str}/modelSegments.called.seg


# Running this function in parallel mode is not allowed
# need python package "gcnvkerbel" installed
..CNV_Germline_Cohort <- function(bam.files.c, gatk, sample.each.time, ref.fasta, ref.dict, interval, blacklist_interval, segmental_duplication_track_bed, mappability_track_bed, contig.ploidy.priors.table) {





   # ========================DetermineGermlineContigPloidy========================
   # r = 3, gcnvkernel is not installed




}



# run in parallel mode
# named.vector = c(R1 = '~/projects/test/fastq/1801169-RT_R1.fq', R2 = '~/projects/test/fastq/1801169-RT_R2.fq')
# return 0
# output:
..Preprocessing.RNA <- function(named.vector, gatk, star, samtools, ref.fasta, StarReferencesFolder, known.sites.VCFs.c, dbsnp, scatter.count) {

   Sys.sleep(runif(1, 0.1, 2))

   fastq.1.str = named.vector['R1']
   fastq.2.str = named.vector['R2']
   gatk = gatk
   ref.fasta = ref.fasta
   StarReferencesFolder = StarReferencesFolder
   star = star
   known.sites.VCFs = known.sites.VCFs.c
   dbsnp = dbsnp
   scatter.count = scatter.count

   sample.name.str = str_extract(basename(fastq.1.str), pattern = '^\\S+(?=_R1)')  # 1801169-RT
   if (sample.name.str != str_extract(basename(fastq.2.str), pattern = '^\\S+(?=_R2)')) {
      message(sprintf('%s and %s have different sample names.', fastq.1.str, fastq.2.str))
      return(1)
   }

   temp.folder = sprintf('temp/preprocess.rna/%s', sample.name.str)
   unlink(temp.folder, recursive = T)   # important
   dir.create(temp.folder, recursive = T, showWarnings = F)   # temp folder, must be empty

   output.folder = 'bamfiles'
   dir.create(output.folder, showWarnings = F)         # output folder

   #===========================STARAlign                   output: {temp.folder}/Aligned.sortedByCoord.out.bam===========================
   if (str_detect(fastq.1.str, pattern = '.gz$')) {
      input.gzip.str = '--readFilesCommand "gzip -cd"'
   } else {
      input.gzip.str = ''
   }

   read.length.int = get_fastq_read_length(fastq.1.str)

   command.str = sprintf('%s --genomeDir %s --runThreadN %s --readFilesIn %s %s %s --sjdbOverhang %s --outFileNamePrefix %s/ --outSAMtype BAM SortedByCoordinate --twopassMode Basic --chimOutType Junctions --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30', star, StarReferencesFolder, scatter.count, fastq.1.str, fastq.2.str, input.gzip.str, read.length.int - 1, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: STARAlign'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/Aligned.sortedByCoord.out.bam
   # output: {temp.folder}/Chimeric.out.junction

   #===========================AddGroupID and Fix tag      output: {temp.folder}/Sec.rem.bam===========================
   command.str = sprintf('%s AddOrReplaceReadGroups -I %s/Aligned.sortedByCoord.out.bam -O %s/Add.out.bam --RGID %s --RGSM %s --RGLB SureSelect --RGPL ILLUMINA --RGPU READ_LENGTH_%s', gatk, temp.folder, temp.folder, sample.name.str, sample.name.str, read.length.int)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: AddOrReplaceReadGroups'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/Add.out.bam

   command.str = sprintf('%s SetNmMdAndUqTags -I %s/Add.out.bam -O %s/Fix.out.bam -R %s', gatk, temp.folder, temp.folder, ref.fasta)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: SetNmMdAndUqTags'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/Fix.out.bam

   command.str = sprintf('%s view -bh -F 0x900 -@ %s -o %s/Sec.rem.bam %s/Fix.out.bam', samtools, scatter.count, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: Remove Secondary'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/Sec.rem.bam

   #===========================Remove duplicates           output: {temp.folder}/Dup.rem.bam===========================
   command.str = sprintf('%s MarkDuplicates -I %s/Sec.rem.bam -O %s/Dup.rem.bam --METRICS_FILE %s/MarkDuplicates.metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT', gatk, temp.folder, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: MarkDuplicates'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/Dup.rem.bam

   #===========================SplitNCigarReads            output: {temp.folder}/Split.cigar.bam===========================
   command.str = sprintf('%s SplitNCigarReads -R %s -I %s/Dup.rem.bam -O %s/Split.cigar.bam', gatk, ref.fasta, temp.folder, temp.folder)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: SplitNCigarReads'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/Split.cigar.bam

   #===========================BaseRecalibrator            output: {temp.folder}/recal_data.table===========================
   input.sites.str = paste0(' --known-sites ', known.sites.VCFs, collapse = '')

   command.str = sprintf('%s --java-options -Xms4g BaseRecalibrator -R %s -I %s/Split.cigar.bam -O %s/recal_data.table --use-original-qualities --known-sites %s %s', gatk, ref.fasta, temp.folder, temp.folder, dbsnp, input.sites.str)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: BaseRecalibrator'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/recal_data.table

   #===========================ApplyBQSR                   output: {temp.folder}/recalibrated.bam===========================
   command.str = sprintf('%s --java-options -Xms3g ApplyBQSR -I %s/Split.cigar.bam --bqsr-recal-file %s/recal_data.table -O %s/recalibrated.bam -R %s --add-output-sam-program-record --use-original-qualities --create-output-bam-index', gatk, temp.folder, temp.folder, temp.folder, ref.fasta)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: ApplyBQSR'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/recalibrated.bam

   #===========================SortSam                     output: {output.folder}/{sample.name.str}.bam===========================
   command.str = sprintf('%s --java-options -Xms8g SortSam --SORT_ORDER coordinate -I %s/recalibrated.bam -O %s/%s.bam --CREATE_INDEX true', gatk, temp.folder, output.folder, sample.name.str)
   message(sprintf('\n[%s] %s: %s', Sys.time(), sample.name.str, 'RNA Preprocess: SortSam'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      unlink(temp.folder, recursive = T)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {output.folder}/{sample.name.str}.bam

   return(0)
}

# run in parallal mode
# bam.interval.c = c(bam.file = 'bamfiles/1234567-RN.bam', interval = 'temp/germline.snp.RNA/temp_NNNN_of_N/scattered.interval_list')
# bam.interval.c, gatk, ref.fasta, dbsnp
..RNAseqGermline.HaplotypeCaller <- function(bam.interval.c, gatk, ref.fasta, dbsnp) {

   bam.file = bam.interval.c['bam.file']
   interval = bam.interval.c['interval']

   prefix = str_extract(interval, pattern = '(?<=temp_)\\d+')

   sample.name.str = file_path_sans_ext(basename(bam.file))

   temp.folder = sprintf('temp/germline.snp.RNA/%s', sample.name.str)
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = sprintf('temp/germline.snp.RNA/%s', sample.name.str)
   dir.create(output.folder, recursive = T, showWarnings = F)

   # ==========================HaplotypeCaller     output: {temp.folder}/RNA.{prefix}.germline.vcf.gz==========================
   command.str = sprintf('%s --java-options -Xms6g HaplotypeCaller -R %s -I %s -L %s -O %s/RNA.%s.germline.vcf.gz --dbsnp %s --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20', gatk, ref.fasta, bam.file, interval, temp.folder, prefix, dbsnp)
   message(sprintf('\n[%s] %s(%s): %s', Sys.time(), sample.name.str, prefix, 'RNA Germline: HaplotypeCaller'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/RNA.{prefix}.germline.vcf.gz

   # ==========================VariantFiltration     output: {temp.folder}/RNA.{prefix}.filtered.germline.vcf.gz==========================
   command.str = sprintf('%s VariantFiltration -R %s -V %s/RNA.%s.germline.vcf.gz -O %s/RNA.%s.filtered.germline.vcf.gz --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0"', gatk, ref.fasta, temp.folder, prefix, output.folder, prefix)
   message(sprintf('\n[%s] %s(%s): %s', Sys.time(), sample.name.str, prefix, 'RNA Germline: VariantFiltration'))

   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      return(0)
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/RNA.{prefix}.filtered.germline.vcf.gz

}
# output: temp/germline.snp.RNA/{sample.name.str}/RNA.{NNNN}.filtered.germline.vcf.gz


# run in parallele mode
# temp.folder = sprintf('temp/somatic.snp/%s', sample.name.str)
# named.vector = c(tumor.bam = '~/projects/test/bamfiles/1801167-PT.bam', normal.bam = '~/projects/test/bamfiles/1801167-PN.bam') normal.bam must not be NA
# strelka.path = '', the command will be <strelka.path>/configManta.py .....
# bed.file = '~/db/BEDs/Agilent_SureSelect_Human_All_Exon_V6_hg19_noXY_noChr.bed.gz'
# ref.fasta = '~/db/mutect2_support/b37/human_g1k_v37.fasta'
#
# temp.folder: "temp/somatic.snp.strelka/1801169"
# output: {output.folder}/{sample.name.str}.somatic.indels.vcf.gz  / .somatic.snvs.vcf.gz
# return 0
..manta.and.strelka.somatic.workflows <- function(named.vector, strelka.path, bcftools.path, vcf2avinput.path, bed.file, ref.fasta, exome, manta.extra.arg = '', strelka.extra.arg = '') {
   library(tools)

   tumor.bam = normalizePath(named.vector['tumor.bam'])
   normal.bam = normalizePath(named.vector['normal.bam'])
   sample.name.str = get_sample_name(tumor.bam)  # 1801169

   if (is.na(tumor.bam)) {
      message('No tumor for this item. Skip')
      return('No tumor for this item')
   }

   if (is.na(normal.bam)) {
      message('No normal for this item. Skip')
      return('No normal for this item')
   }

   if (strelka.path == '') {
      configManta = 'configManta.py'
      configStrelka = 'configureStrelkaSomaticWorkflow.py'
   } else {
      configManta = sprintf('%s/configManta.py', strelka.path)
      configStrelka = sprintf('%s/configureStrelkaSomaticWorkflow.py', strelka.path)
   }

   if (bcftools.path == '') {
      bcftools = 'bcftools'
   } else {
      bcftools = sprintf('%s/bcftools', bcftools.path)
   }

   if (vcf2avinput.path == '') {
      vcf2avinput = 'vcf2avinput.py'
   } else {
      vcf2avinput = sprintf('%s/vcf2avinput.py', vcf2avinput.path)
   }

   ref.fasta = ref.fasta

   temp.folder = sprintf('temp/somatic.snp.strelka/%s', sample.name.str) # "temp/somatic.snp.strelka/1801169"
   unlink(temp.folder, force = T, recursive = T)
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = './somatic.snp.strelka'
   dir.create(output.folder, recursive = T, showWarnings = F)



   Sys.sleep(runif(1, 0.5, 6))
   #message('sleeping...')
   # ============Manta          output: {temp.folder}/manta_workflows/results/variants/candidateSmallIndels.vcf.gz (tbi)==========================
   command.str = sprintf('%s --tumorBam %s --normalBam %s --referenceFasta %s --callRegions %s --runDir %s/manta_workflows %s %s', configManta, tumor.bam, normal.bam, ref.fasta, bed.file, temp.folder, ifelse(exome, yes = '--exome', no = ''), manta.extra.arg)
   r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Config Manta', verbose = F)
   if (r != 0) {return(r)}

   #output: {temp.folder}/manta_workflows/

   command.str = sprintf('%s/manta_workflows/runWorkflow.py', temp.folder)
   r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Run Manta')
   if (r != 0) {return(r)}

   #output: {temp.folder}/manta_workflows/results/variants/candidateSmallIndels.vcf.gz (tbi)


   # ============Strelka        output: {temp.folder}/strelka_workflows/results/variants/somatic.indels.vcf.gz (tbi) / somatic.snvs.vcf.gz==========================
   command.str = sprintf('%s --tumorBam %s --normalBam %s --referenceFasta %s --indelCandidates %s/manta_workflows/results/variants/candidateSmallIndels.vcf.gz --callRegions %s --runDir %s/strelka_workflows %s', configStrelka, tumor.bam, normal.bam, ref.fasta, temp.folder, bed.file, temp.folder, ifelse(exome, yes = '--exome', no = ''), strelka.extra.arg)
   r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Config Strelka', verbose = F)
   if (r != 0) {return(r)}
   #output: {temp.folder}/strelka_workflows/

   command.str = sprintf('%s/strelka_workflows/runWorkflow.py -m local', temp.folder)
   r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Run Strelka')
   if (r != 0) {return(r)}
   # output: {temp.folder}/strelka_workflows/results/variants/somatic.indels.vcf.gz  / somatic.snvs.vcf.gz

   # ============Copy results   output: {output.folder}/{sample.name.str}.somatic.indels.vcf.gz  / .somatic.snvs.vcf.gz==========================
   command.str = sprintf('cp %s/strelka_workflows/results/variants/somatic.indels.vcf.gz %s/%s.somatic.indels.vcf.gz', temp.folder, output.folder, sample.name.str)
   r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Copy results', verbose = F)
   if (r != 0) {return(r)}

   # command.str = sprintf('cp %s/strelka_workflows/results/variants/somatic.indels.vcf.gz.tbi %s/%s.somatic.indels.vcf.gz.tbi', temp.folder, output.folder, sample.name.str)
   # r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Copy results', verbose = F)
   # if (r != 0) {return(r)}

   command.str = sprintf('cp %s/strelka_workflows/results/variants/somatic.snvs.vcf.gz %s/%s.somatic.snvs.vcf.gz', temp.folder, output.folder, sample.name.str)
   r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Copy results', verbose = F)
   if (r != 0) {return(r)}

   # command.str = sprintf('cp %s/strelka_workflows/results/variants/somatic.snvs.vcf.gz.tbi %s/%s.somatic.snvs.vcf.gz.tbi', temp.folder, output.folder, sample.name.str)
   # r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'Copy results', verbose = F)
   # if (r != 0) {return(r)}
   # output: {output.folder}/{sample.name.str}.somatic.indels.vcf.gz / .somatic.snvs.vcf.gz

   # ============Normalize the result   Output: {temp.folder}/strelka_workflows/results/variants/somatic.indels.vcf.gz (tbi) / somatic.snvs.vcf.gz==========================
   r = system(sprintf('%s --version', bcftools), ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      # ~/bin/bcftools-1.11/bcftools norm --force -m-both -o ex1.step1.vcf ex1.vcf.gz
      input = sprintf('%s/strelka_workflows/results/variants/somatic.indels.vcf.gz', temp.folder)
      output = sprintf('%s/strelka_workflows/results/variants/temp.vcf.gz', temp.folder)
      command.str = sprintf('%s norm --force -m-both -o %s %s', bcftools, output, input)
      r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'bcftools split indels', verbose = T)
      if (r != 0) {return(r)}

      # ~/bin/bcftools-1.11/bcftools norm --force -f ~/db/mutect2_support/b37/human_g1k_v37.fasta -o ex1.step2.vcf ex1.step1.vcf
      input = output
      output = sprintf('%s/strelka_workflows/results/variants/normalized.indels.vcf', temp.folder)
      command.str = sprintf('%s norm --force -f %s -o %s %s', bcftools, ref.fasta, output, input)
      r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'bcftools norm indels', verbose = T)
      if (r != 0) {return(r)}

      # ~/bin/vcf2avinput.py -s "TUMOR" -o normalized.somatic.snvs.avinput normalized.somatic.snvs.vcf
      input = output
      output = sprintf('%s/strelka_workflows/results/variants/normalized.indels.avinput', temp.folder)
      command.str = sprintf('%s -s "TUMOR" -o %s %s', vcf2avinput, output, input)
      r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'convert avinput indels', verbose = T)
      if (r != 0) {return(r)}

      indel.df = read.csv(output, sep = '\t', header = F, stringsAsFactors = F)
      indel.df$V6 = sample.name.str

      # # # # # # # # # # # # # # # # # #
      # ~/bin/bcftools-1.11/bcftools norm --force -m-both -o ex1.step1.vcf ex1.vcf.gz
      input = sprintf('%s/strelka_workflows/results/variants/somatic.snvs.vcf.gz', temp.folder)
      output = sprintf('%s/strelka_workflows/results/variants/temp.vcf.gz', temp.folder)
      command.str = sprintf('%s norm --force -m-both -o %s %s', bcftools, output, input)
      r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'bcftools split snvs', verbose = T)
      if (r != 0) {return(r)}

      # ~/bin/bcftools-1.11/bcftools norm --force -f ~/db/mutect2_support/b37/human_g1k_v37.fasta -o ex1.step2.vcf ex1.step1.vcf
      input = output
      output = sprintf('%s/strelka_workflows/results/variants/normalized.snvs.vcf', temp.folder)
      command.str = sprintf('%s norm --force -f %s -o %s %s', bcftools, ref.fasta, output, input)
      r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'bcftools norm snvs', verbose = T)
      if (r != 0) {return(r)}

      # ~/bin/vcf2avinput.py -s "TUMOR" -o normalized.somatic.snvs.avinput normalized.somatic.snvs.vcf
      input = output
      output = sprintf('%s/strelka_workflows/results/variants/normalized.snvs.avinput', temp.folder)
      command.str = sprintf('%s -s "TUMOR" -o %s %s', vcf2avinput, output, input)
      r = run.single.command(command = command.str, sample.name = sample.name.str, abbr = 'convert avinput snvs', verbose = T)
      if (r != 0) {return(r)}

      snv.df = read.csv(output, sep = '\t', header = F, stringsAsFactors = F)
      snv.df$V6 = sample.name.str

      total.df = rbind(snv.df, indel.df)
      total.df = total.df[order(total.df$V1, total.df$V2, total.df$V3, total.df$V4, total.df$V5), ]

      output = sprintf('%s/%s.avinput', output.folder, sample.name.str)
      write.table(total.df, output, sep = '\t', quote = F, row.names = F, col.names = F)
   }

   unlink(temp.folder, recursive = T)
   return(0)

}



# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# fastq.list is a text file. Each line contains a fastq file in absolute path.
# The fastq files pair must named as <sample_name>_R1.<ext>, <sample_name>_R2.<ext>, the script use '_R1.' and '_R2.' to recognize sample name and pair
# The parameter.file is a headerless tab-delimited text file with two columns. First column is the parameter name and the second is content. Check the example file to see detail
# Return processed bam file names in a vector and a bam file list
GATK.Preprocess <- function(fastq.list, parameter.file) {

   # read parameters
   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   mapper = parameters.df['mapper', 'V2']
   mapper.index = parameters.df['mapper.index', 'V2']
   samtools = parameters.df['samtools', 'V2']
   gatk = parameters.df['gatk', 'V2']
   scatter.count = as.integer(parameters.df['scatter.count', 'V2'])
   sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])
   keep.unBQSR = as.logical(parameters.df['keep.unBQSR', 'V2'])

   ref.fasta = parameters.df['ref.fasta', 'V2']
   ref.dict = parameters.df['ref.dict', 'V2']
   known.sites.VCFs.c = str_split(parameters.df['known.sites.VCFs', 'V2'], pattern = '\\|', simplify = T)[1, ]


   output.folder = 'bamfiles'

   # read fastq.list
   fastq.files.list = generate.fastq.list(fastq.list)

   # named.vector = c(R1 = '~/projects/test/fastq/1801169-DT_R1.fq', R2 = '~/projects/test/fastq/1801169-DT_R2.fq')
   r = parallel::mclapply(fastq.files.list, FUN = ..Preprocessing.DNA, mapper = mapper, mapper.index = mapper.index, samtools = samtools, gatk = gatk, scatter.count = scatter.count, keep.unBQSR = keep.unBQSR, ref.fasta = ref.fasta, ref.dict = ref.dict, known.sites.VCFs = known.sites.VCFs.c, mc.cores = sample.each.time)
   print(r)

   # generate bam list
   files.c = sprintf('%s/%s.bam', output.folder, names(fastq.files.list))
   message(sprintf('Generate %s bam files.', sum(file.exists(files.c))))
   if (sum(file.exists(files.c)) != length(fastq.files.list)) {
      message('Bam files count is not equal to sample names count, something wrong ?')
   }

   files.c = normalizePath(files.c[file.exists(files.c)])
   write.table(data.frame(V1 = files.c), sprintf('%s/bams.DNA.tsv', output.folder), sep = '\t', row.names = F, col.names = F, quote = F)
   return(0)
}


# bams.list is a text file. Each line contains a bam file in absolute path.
# The bam files pair must named as <sample_name>-DT.bam, <sample_name>-DN.bam, the script use '-DT.' ('-RT.') and '-DN.' ('-RN.') to recognize sample name and tumor/normal
# The parameter.file is a headerless tab-delimited text file with two columns. First column is the parameter name and the second is content. Check the example file to see detail
# Return processed vcf file names in a vector and creates a vcf file list
GATK.Somatic.Variants <- function(bams.list, parameter.file) {

   # read parameters
   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   gatk = parameters.df['gatk', 'V2']
   gatk.safe = parameters.df['gatk.safe', 'V2']
   scatter.count = as.integer(parameters.df['scatter.count', 'V2'])
   sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])
   interval = parameters.df['interval.file', 'V2']
   ref.fasta = parameters.df['ref.fasta', 'V2']
   ref.dict = parameters.df['ref.dict', 'V2']
   variants.for.contamination = parameters.df['variants.for.contamination', 'V2']

   pon = parameters.df['pon', 'V2']
   gnomad = parameters.df['gnomad', 'V2']
   realignment.index.bundle = parameters.df['realignment.index.bundle', 'V2']
   gga.vcf = parameters.df['gga.vcf', 'V2']
   m2.extra.args = parameters.df['m2.extra.args', 'V2']
   getpileupsummaries.extra.args = parameters.df['getpileupsummaries.extra.args', 'V2']
   filter.mutect2.extra.args = parameters.df['filter.mutect2.extra.args', 'V2']
   realignment.extra.args = parameters.df['realignment.extra.args', 'V2']

   # read bams.list
   sample.list = generate.bam.list(bams.list)

   temp.folder = 'temp/somatic.snp'
   dir.create(temp.folder, recursive = T, showWarnings = F) # "temp/somatic.snp/1801169-PT/"

   output.folder = 'somatic.snp'
   dir.create(output.folder, recursive = T, showWarnings = F)

   # ==========================SplitIntervals  output: {temp.folder}/0000-scattered.interval_list .... 0007-scattered.interval_list==========================
   interval.files.c = list.files(temp.folder, pattern = '.interval_list$', full.names = T)
   unlink(interval.files.c, force = T)

   command.str = sprintf('%s SplitIntervals -R %s -L %s -scatter %s -O %s', gatk, ref.fasta, interval, scatter.count, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'Somatic SplitIntervals'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit(save = 'no', status = r)
   }
   interval.files.c = list.files(temp.folder, pattern = '.interval_list$', full.names = T)
   # output:  {temp.folder}/0000-scattered.interval_list .... 0007-scattered.interval_list

   #===========================Index PON       output: pon.vcf.idx =========================
   if (pon != '') {
      command.str = sprintf('%s IndexFeatureFile -I %s', gatk, pon)
      message(sprintf('\n[%s] %s', Sys.time(), 'IndexFeatureFile PON'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }
   }

   # ==========================Mutect2         output: somatic.snp/<sample_name>.vcf==========================
   r = parallel::mclapply(sample.list, FUN = ..Mutect2, gatk = gatk, gatk.safe = gatk.safe, intervals = interval.files.c, scatter.count = scatter.count, ref.fasta = ref.fasta, ref.dict = ref.dict, variants.for.contamination = variants.for.contamination, pon = pon, gnomad = gnomad, gga.vcf = gga.vcf, realignment.index.bundle = realignment.index.bundle, m2.extra.args = m2.extra.args, getpileupsummaries.extra.args = getpileupsummaries.extra.args, filter.mutect2.extra.args = filter.mutect2.extra.args, realignment.extra.args = realignment.extra.args, mc.cores = sample.each.time)
   message()
   print(r)
   message(r)
}

# output: {output.folder}/somatic.pon.vcf
GATK.Somatic.PON <- function(bams.list, parameter.file) {

   # read parameters
   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   gatk = parameters.df['gatk', 'V2']
   gatk.safe = parameters.df['gatk.safe', 'V2']
   scatter.count = as.integer(parameters.df['scatter.count', 'V2'])
   sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])
   interval = parameters.df['interval.file', 'V2']
   ref.fasta = parameters.df['ref.fasta', 'V2']
   ref.dict = parameters.df['ref.dict', 'V2']
   variants.for.contamination = parameters.df['variants.for.contamination', 'V2']

   gnomad = parameters.df['gnomad', 'V2']
   realignment.index.bundle = parameters.df['realignment.index.bundle', 'V2']
   gga.vcf = parameters.df['gga.vcf', 'V2']
   getpileupsummaries.extra.args = parameters.df['getpileupsummaries.extra.args', 'V2']
   filter.mutect2.extra.args = parameters.df['filter.mutect2.extra.args', 'V2']
   realignment.extra.args = parameters.df['realignment.extra.args', 'V2']
   create_pon_extra_args = parameters.df['create_pon_extra_args', 'V2']

   mutect2 = as.logical(parameters.df['mutect2', 'V2'])
   create.panel = as.logical(parameters.df['create.panel', 'V2'])

   temp.folder = 'temp/somatic.snp'
   message(sprintf('Deleting %s', temp.folder))
   dir.create(temp.folder, recursive = T, showWarnings = F) # "temp/somatic.snp"

   output.folder = 'somatic.snp'  # the normal samples vcf filse must be here
   dir.create(output.folder, recursive = T, showWarnings = F)

   # the same as ..Mutect2, except the following:
   # add: bams.list, create_pon_extra_args
   # gnomad, optional -> positional
   # removed: pon, m2.extra.args
   # return(0)

   r = system('ulimit -c unlimited', ignore.stdout = T, ignore.stderr = T)

   # ==========================SplitIntervals  output: {temp.folder}/0000-scattered.interval_list .... NNNN-scattered.interval_list==========================
   interval.files.c = list.files(temp.folder, pattern = '.interval_list$', full.names = T)
   unlink(interval.files.c)

   command.str = sprintf('%s SplitIntervals -R %s -L %s -scatter %s -O %s', gatk, ref.fasta, interval, scatter.count, temp.folder)
   message(sprintf('\n[%s] %s', Sys.time(), 'Somatic PON SplitIntervals'))
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      quit(save = 'no', status = r)
   }
   # output: {temp.folder}/0000-scattered.interval_list .... NNNN-scattered.interval_list

   # ==========================Mutect2         output: {output.folder}/<sample_name>-DN.vcf.gz==========================
   #> sample.list
   #$`1801167`
   #tumor.bam
   #"~/projects/test/bamfiles/1801167-PN.bam"
   #$`1801169`
   #tumor.bam
   #"~/projects/test/bamfiles/1801169-PN.bam"

   # import all normal samples into a list, discard the tumor samples
   # We have to change the vector name from 'normal.bam' to 'tumor.bam'
   sample.list = generate.bam.list(bams.list, allow.type = 'normal')
   for (i in 1:length(sample.list)) {
      names(sample.list[[i]]) = 'tumor.bam'
   }

   interval.files.c = list.files(temp.folder, pattern = '.interval_list$', full.names = T)

   if (mutect2) {
      r = parallel::mclapply(sample.list, FUN = ..Mutect2, gatk = gatk, gatk.safe = gatk.safe, intervals = interval.files.c, scatter.count = scatter.count, ref.fasta = ref.fasta, ref.dict = ref.dict, variants.for.contamination = variants.for.contamination, gnomad = gnomad, gga.vcf = gga.vcf, realignment.index.bundle = realignment.index.bundle, m2.extra.args = '--max-mnp-distance 0', getpileupsummaries.extra.args = getpileupsummaries.extra.args, filter.mutect2.extra.args = filter.mutect2.extra.args, realignment.extra.args = realignment.extra.args, mc.cores = sample.each.time)

      if (all(r == 0)) {
         message(sprintf('\n[%s] SUCCESS: Mutect2', Sys.time()))
      } else {
         print(r)
         message('Mutect2_Panel_M2: ', r)
         quit(save = 'no', status = 1)
      }
   }
   # output: {output.folder}/<sample_name>-DN.vcf.gz


   #===========================CreatePanel     output: {output.folder}/somatic.pon.vcf=======================
   ..CreatePanel <- function(interval, vcf.files.c, create_pon_extra_args) {
      Sys.sleep(runif(1, 0.1, 2))

      prefix = str_extract(basename(interval), pattern = '\\d+')   # NNNN
      message(sprintf('\n[%s] %s: %s', Sys.time(), prefix, 'Deleting existed genomicsdb folder ...'))
      unlink(sprintf('%s/genomicsdb_%s', temp.folder, prefix), recursive = T, force = T)

      # GenomicsDBImport
      dir.create(sprintf('%s/GenomicsDBImport.temp/%s', temp.folder, prefix), recursive = T, showWarnings = F) # "{temp.folder}/GenomicsDBImport.temp/NNNN"

      command.str = sprintf('%s --java-options -Xmx32g GenomicsDBImport --tmp-dir %s/GenomicsDBImport.temp/%s --genomicsdb-workspace-path %s/genomicsdb_%s -R %s -L %s %s %s --genomicsdb-shared-posixfs-optimizations', gatk, temp.folder, prefix, temp.folder, prefix, ref.fasta, interval, paste0(' -V ', vcf.files.c, collapse = ''), create_pon_extra_args)
      message(sprintf('\n[%s] %s: %s', Sys.time(), prefix, 'GenomicsDBImport'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }

      unlink(sprintf('%s/GenomicsDBImport.temp/%s', temp.folder, prefix))   # delete GenomicsDBImport temp folder
      # output: {temp.folder}/genomicsdb_NNNN

      # CreateSomaticPanelOfNormals
      command.str = sprintf('%s --java-options -Xmx32g CreateSomaticPanelOfNormals -R %s -V gendb://%s/genomicsdb_%s -O %s/pon_%s_raw.vcf --germline-resource %s %s', gatk, ref.fasta, temp.folder, prefix, temp.folder, prefix, gnomad, create_pon_extra_args)
      message(sprintf('\n[%s] %s: %s', Sys.time(), prefix, 'CreateSomaticPanelOfNormals'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }

      unlink(sprintf('%s/genomicsdb_%s', temp.folder, prefix), recursive = T, force = T)
      #output: {temp.folder}/pon_NNNN_raw.vcf

      # SortVcf, sometime the pon_NNNN vcf need to be sorted
      command.str = sprintf('%s --java-options -Xmx4g SortVcf -I %s/pon_%s_raw.vcf -O %s/pon_%s.vcf', gatk, temp.folder, prefix, temp.folder, prefix)
      message(sprintf('\n[%s] %s: %s', Sys.time(), prefix, 'SortVcf'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit(save = 'no', status = r)
      }
      #output: {temp.folder}/pon_NNNN.vcf

   }
   #output: {temp.folder}/pon_NNNN.vcf

   interval.files.c = list.files(temp.folder, pattern = '.interval_list$', full.names = T)

   if (create.panel) {

      # generate vcf files vector
      bam.files.c = generate.bam.list(bams.list, allow.type = 'normal', return.type = 'vector')
      vcf.files.c = c()
      for (i in 1:length(bam.files.c)) {
         temp = tools::file_path_sans_ext(basename(bam.files.c[i]))  # 1801167-XN
         vcf.file = sprintf('%s/%s.vcf.gz', output.folder, temp)
         message(sprintf('Try to import %s (%s/%s)', vcf.file[1], i, length(bam.files.c)))
         if (file.exists(vcf.file)) {
            vcf.files.c = c(vcf.files.c, vcf.file)
            message(sprintf('Import %s', vcf.file))
         }
      }
      message(sprintf('%s VCF files.', length(vcf.files.c)))

      # create scattered panel vcf files
      r = parallel::mclapply(interval.files.c, FUN = ..CreatePanel, vcf.files.c = vcf.files.c, create_pon_extra_args = '', mc.cores = min(sample.each.time, floor(get_machine_mem()/32)))
      message()
      message('CreatePanel:')
      print(paste0(interval.files.c, ' (', r, ')'))
      message(r)

      if (all(r == 0)) {
         message(sprintf('\n[%s] SUCCESS: CreatePanel', Sys.time()))
      } else {
         print(r)
         message('CreatePanel: ', r)
         quit(save = 'no', status = 1)
      }
      # output: {temp.folder}/pon_0000.vcf, .... {temp.folder}/pon_NNNN.vcf

      # mergeVCFs
      pon.vcf.c = paste0(temp.folder, '/pon_', str_extract(basename(interval.files.c), pattern = '\\d+'), '.vcf')  # {temp.folder}/pon_NNNN.vcf
      input.c = paste0(' -I ', pon.vcf.c, collapse = '')
      command.str = sprintf('%s --java-options -Xmx8g MergeVcfs %s -O %s/pon.somatic.vcf', gatk, input.c, temp.folder)
      message(sprintf('\n[%s] %s', Sys.time(), 'MergeVcfs'))
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         file.copy(sprintf('%s/pon.somatic.vcf', temp.folder), sprintf('%s/somatic.pon.vcf', output.folder))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit('no', status = r)
      }
      # output: {output.folder}/somatic.pon.vcf

      if (r == 0) {
         message('Panel vcf file created.')
         return(r)
      }

   }

}

GATK.Germline.Variants <- function(bams.list, parameter.file) {

      # read parameters
      parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

      haplotypeCaller = as.logical(parameters.df['haplotypeCaller', 'V2'])
      joint.calling = as.logical(parameters.df['joint.calling', 'V2'])
      gatk = parameters.df['gatk', 'V2']
      gatk.safe = parameters.df['gatk.safe', 'V2']
      scatter.count = as.integer(parameters.df['scatter.count', 'V2'])
      sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])
      interval = parameters.df['unpadded_intervals_file', 'V2']
      ref.fasta = parameters.df['ref.fasta', 'V2']

      max.gaussians.indel = as.integer(parameters.df['max.gaussians.indel', 'V2'])
      max.gaussians.snp = as.integer(parameters.df['max.gaussians.snp', 'V2'])
      str.table = parameters.df['str.table', 'V2']
      dbsnp = parameters.df['dbsnp', 'V2']
      mills.resource = parameters.df['mills.resource', 'V2']
      axiomPoly.resource = parameters.df['axiomPoly.resource', 'V2']
      hapmap.resource = parameters.df['hapmap.resource', 'V2']
      omni.resource = parameters.df['omni.resource', 'V2']
      one.thousand.genomes.resource = parameters.df['one.thousand.genomes.resource', 'V2']

      bam.files.c = generate.bam.list(bams.list, return.type = 'vector')

      # ================================haplotypeCaller       output: {temp.folder}/{sample.name}/HaplotypeCaller.g.vcf.gz .......===============================
      # catalogue the tumor samples and normal samples into list
      # ..HaplotypeCallerGvcf_GATK4 <- function(bam.file, gatk, interval, ref.fasta)
      if (haplotypeCaller) {
         r = parallel::mclapply(bam.files.c, FUN = ..HaplotypeCallerGvcf, gatk = gatk, interval = interval, ref.fasta = ref.fasta, str.table = str.table, mc.cores = round(get_cpu_count()/3))
         message('\nHaplotypeCallerGvcf reuslt: \n')
         print(paste0(bam.files.c, ' (', r, ')'))
      }
      # output: {temp.folder}/{sample.name}/HaplotypeCaller.g.vcf.gz .......

      # ================================Jointcalling          output: {output.folder}/jointGenotyping.vcf===============================
      # temp/variant_germline/1801169-PN/HaplotypeCaller.g.vcf.gz
      # ..JointGenotyping <- function(gVCFs.files.c, gatk, interval, scatter.count, sample.each.time, ref.fasta, dbsnp, mills.resource, axiomPoly.resource, hapmap.resource, omni.resource, one.thousand.genomes.resource, indel.max.gaussians = 4, snp.max.gaussians = 6)
      gVCFs.files.c = c()
      for (sample.name.str in file_path_sans_ext(basename(bam.files.c))) {
         gVCFs.file = sprintf('temp/germline.snp/%s/HaplotypeCaller.g.vcf.gz', sample.name.str)
         if (file.exists(gVCFs.file) & file.size(gVCFs.file) != 0) {
            gVCFs.files.c = c(gVCFs.files.c, gVCFs.file)
         } else {
            message(sprintf('Can not find %s, skip.', gVCFs.file))
         }
      }

      if (length(gVCFs.files.c) == 0) {
         message('No gVCF files found')
         return(1)
      } else {
         message('Found following gVCF files:')
         print(gVCFs.files.c)
         message(sprintf('%s gVCF files.', length(gVCFs.files.c)))
      }

      if (joint.calling) {
         r = ..JointGenotyping(gVCFs.files.c, gatk, gatk.safe, interval, scatter.count, sample.each.time, ref.fasta, dbsnp, mills.resource, axiomPoly.resource, hapmap.resource, omni.resource, one.thousand.genomes.resource, max.gaussians.indel, max.gaussians.snp)
         message('\nJoint Genotyping: ', r)
         quit(save = 'no', status = r)
      }
      # output: {output.folder}/jointGenotyping.vcf
}

# Generate segment file for GISTIC2
# mkdir gistic2_output_GATK
# gistic2 -b gistic2_output_GATK -seg segmentation_file_for_gistic2.txt -refgene ~/db/GISTIC2.hg19.mat -conf 0.9
GATK.Somatic.CNV <- function(bams.list, parameter.file) {

   # read parameters
   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   gatk = parameters.df['gatk', 'V2']
   sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])
   ref.fasta = parameters.df['ref.fasta', 'V2']
   ref.dict = parameters.df['ref.dict', 'V2']
   interval = parameters.df['interval.file', 'V2']
   create.panel = as.logical(parameters.df['create.panel', 'V2'])
   call.somatic.cnv = as.logical(parameters.df['call.somatic.cnv', 'V2'])

   common.sites = parameters.df['common.sites', 'V2']
   blacklist_interval = parameters.df['blacklist_interval', 'V2']
   segmental_duplication_track_bed = parameters.df['segmental_duplication_track_bed', 'V2']
   mappability_track_bed = parameters.df['mappability_track_bed', 'V2']

   gistic2 = parameters.df['gistic2', 'V2']
   gistic2_refgene_mat = parameters.df['gistic2_refgene_mat', 'V2']
   gistic2_conf = as.numeric(parameters.df['gistic2_conf', 'V2'])

   temp.folder = 'temp/somatic.cnv'
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = 'somatic.CNV'
   dir.create(output.folder, recursive = T, showWarnings = F)

   # =========================generate CNV PON      output: {temp.folder}/cnv.somatic.pon.hdf5===================================
   if (create.panel) {
      message('============Create CNV PON============')
      sample.list = generate.bam.list(bams.list, allow.type = 'normal')

      r = ..CNV_Somatic_Panel(sample.list, gatk, sample.each.time, ref.fasta, ref.dict, interval, common.sites, blacklist_interval, segmental_duplication_track_bed, mappability_track_bed)
      if (r == 0) {
         message(sprintf('\n[%s] CNVSomaticPanelWorkflow SUCCESS', Sys.time()))
      } else {
         message(sprintf('\n[%s] CNVSomaticPanelWorkflow Failed', Sys.time()))
         quit('no', status = r)
      }
      # output: temp/cnv_somatic/cnv.somatic.pon.hdf5
   }


   # =========================call somatic CNV      output: {temp.folder}/{sample.name.str}/modelSegments.called.seg=========================
   if (call.somatic.cnv) {
      preprocessed.interval = sprintf('%s/preprocessed.interval_list', temp.folder)
      if (!file.exists(preprocessed.interval)) {
         message(sprintf('%s not found. Exit', preprocessed.interval))
         quit(save = 'no', status = 1)
      }

      sample.list = generate.bam.list(bams.list, allow.type = 'both')
      r = parallel::mclapply(sample.list, FUN = ..CNV_Somatic_Pair_Call, gatk = gatk, ref.fasta = ref.fasta, common.sites = common.sites, preprocessed.interval = preprocessed.interval, mc.cores = sample.each.time)
      message('\nCNVSomaticPairWorkflow:\n' )
      print(r)
      message(r)
   }

   # =========================merge segments      output: {temp.folder}/segmentation_file_for_gistic2.txt=========================
   # process the calling results and generate the GISTIC2 inputs
   # temp/cnv_somatic/%s/modelSegments.called.seg
   if (!create.panel & !call.somatic.cnv) {quit(save = 'no', status = 0)}

   segment.df = data.frame(matrix(nrow = 0, ncol = 6))
   colnames(segment.df) = c('SAMPLE', "CONTIG", "START", "END", "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO")
   for (i in 1:length(sample.list)) {
      if (is.na(sample.list[[i]]['tumor.bam'])) {next} # if this patient doesn't have tumor sample
      segment.file.str = sprintf('%s/%s/modelSegments.called.seg', temp.folder, file_path_sans_ext(basename(sample.list[[i]]['tumor.bam'])))
      if (!file.exists(segment.file.str)) {
         message(segment.file.str, ' does not exist, skip!')
         next
      }
      message('Import ', segment.file.str)
      temp.df = read.table(file = segment.file.str, row.names = NULL, header = TRUE, sep = "\t", comment.char = "@")
      temp.df = temp.df[, c("CONTIG", "START", "END", "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO")]
      temp.df = data.frame(SAMPLE = file_path_sans_ext(basename(sample.list[[i]]['tumor.bam'])) , temp.df)
      segment.df = rbind(segment.df, temp.df)
   }
   write.table(segment.df, file = sprintf('%s/segmentation_file_for_gistic2.txt', temp.folder), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
   message()
   message(sprintf('Successfully write %s segments to %s/segmentation_file_for_gistic2.txt', i, temp.folder))

   # =========================GISTIC2               output: {output.folder}=========================
   # gistic2 -b gistic2_output_GATK -seg segmentation_file_for_gistic2.txt -refgene ~/db/GISTIC2.hg19.mat -conf 0.9
   command.str = sprintf('%s -b %s -seg %s/segmentation_file_for_gistic2.txt -refgene %s -conf %s', gistic2, output.folder, temp.folder, gistic2_refgene_mat, gistic2_conf)
   if (nrow(segment.df) > 0) {
      dir.create(output.folder, showWarnings = F)
      message(sprintf('\n[%s] %s', Sys.time(), 'GISTIC2'))
      r = system(command.str, ignore.stdout = T, ignore.stderr = F)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         file.rename(sprintf('%s/segmentation_file_for_gistic2.txt', temp.folder), sprintf('%s/segmentation_file_for_gistic2.txt', output.folder))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         quit('no', status = r)
      }
   }

}


# fastq.list is a text file. Each line contains a fastq file in absolute path.
# The fastq files pair must named as <sample_name>_R1.<ext>, <sample_name>_R2.<ext>, the script use '_R1.' and '_R2.' to recognize sample name and pair
GATK.Preprocess.RNA <- function(fastq.list, parameter.file) {

   # read parameters
   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   gatk = parameters.df['gatk', 'V2']
   star = parameters.df['star', 'V2']
   StarReferencesFolder = parameters.df['StarReferencesFolder', 'V2']
   samtools = parameters.df['samtools', 'V2']

   ref.fasta = parameters.df['ref.fasta', 'V2']
   ref.gtf = parameters.df['ref.gtf', 'V2']
   known.sites.VCFs.c = str_split(parameters.df['known.sites.VCFs', 'V2'], pattern = '\\|', simplify = T)[1, ]
   dbsnp = parameters.df['dbsnp', 'V2']
   scatter.count = as.integer(parameters.df['scatter.count', 'V2'])
   sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])

   # read fastq.list
   fastq.files.list = generate.fastq.list(fastq.list)

   read.length.int = get_fastq_read_length(fastq.files.list[[1]]['R1'])

   r = system('ulimit -n unlimited', ignore.stdout = T, ignore.stderr = T)

   output.folder = 'bamfiles'
   # ==========================generate STAR reference if not avaliable==========================
   if (is.na(StarReferencesFolder) | is.null(StarReferencesFolder) | StarReferencesFolder == '') {
      dir.create(StarReferencesFolder, showWarnings = F, recursive = T, mode = '0755')
      command.str = sprintf('%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang %s --runThreadN %s', star, StarReferencesFolder, ref.fasta, ref.gtf, read.length.int - 1, scatter.count*sample.each.time)
      message(sprintf('\n[%s] %s', Sys.time(), 'RNA Preprocess: GenerateReference'))
      dir.create(StarReferencesFolder, showWarnings = F, recursive = T, mode = '0755')
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      }
   }

   # ==========================Preprocess==========================
   # ..Preprocessing.RNA <- function(named.vector, gatk, star, samtools, ref.fasta, StarReferencesFolder, known.sites.VCFs.c, dbsnp, scatter.count)
   r = parallel::mclapply(fastq.files.list, FUN = ..Preprocessing.RNA, gatk = gatk, star = star, samtools = samtools, ref.fasta = ref.fasta, StarReferencesFolder = StarReferencesFolder, known.sites.VCFs.c = known.sites.VCFs.c, dbsnp = dbsnp, scatter.count = scatter.count, mc.cores = sample.each.time)
   message('\nRNAseqGermline.Preprocess:\n' )
   print(r)

   # ==========================generate bam list==========================
   files.c = sprintf('%s/%s.bam', output.folder, names(fastq.files.list))
   message(sprintf('Generate %s bam files.', sum(file.exists(files.c))))
   if (sum(file.exists(files.c)) != length(fastq.files.list)) {
      message('Bam files count is not equal to sample names count, something wrong ?')
   }

   files.c = normalizePath(files.c[file.exists(files.c)])
   write.table(data.frame(V1 = files.c), sprintf('%s/bams.RNA.tsv', output.folder), sep = '\t', row.names = F, col.names = F, quote = F)

}


# gatk
# interval
# scatter.count
GATK.Germline.Variants.RNA <- function(bams.list, parameter.file) {
   # read parameters
   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   gatk = parameters.df['gatk', 'V2']
   interval = parameters.df['interval', 'V2']
   sample.each.time = as.integer(parameters.df['sample.each.time', 'V2'])
   scatter.count = as.integer(parameters.df['scatter.count', 'V2'])

   ref.fasta = parameters.df['ref.fasta', 'V2']
   dbsnp = parameters.df['dbsnp', 'V2']

   # read bam files
   bam.files.c = generate.bam.list(bams.list, return.type = 'vector')

   temp.folder = 'temp/germline.snp.RNA'
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = 'germline.snp.RNA'
   dir.create(output.folder, recursive = T, showWarnings = F)

   # ==========================split interval                          output: {temp.folder}/temp_NNNN_of_{scatter.count}/scattered.interval_list==========================
   command.str = sprintf('%s --java-options -Xms1g IntervalListTools -I %s --SCATTER_COUNT %s -O %s --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE true --SORT true', gatk, interval, scatter.count, temp.folder)
   message('\nRNA Germline: IntervalListTools')
   r = system(command.str, ignore.stderr = T, ignore.stdout = T)
   if (r == 0) {
      message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
   } else {
      message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
      return(1)
   }
   # output: {temp.folder}/temp_NNNN_of_{scatter.count}/scattered.interval_list

   intervals.c = list.files(path = temp.folder, pattern = 'scattered.interval_list', recursive = T, full.names = T)

   # ==========================RNAseqGermline.HaplotypeCaller          output: temp/germline.snp.RNA/{sample.name.str}/RNA.{NNNN}.filtered.germline.vcf.gz==========================
   #..RNAseqGermline.HaplotypeCaller <- function(bam.interval.c, gatk, ref.fasta, dbsnp)
   # bam.interval.c = c(bam.file = 'bamfiles/1234567-RN.bam', interval = 'temp/germline.snp.RNA/temp_NNNN_of_N/scattered.interval_list')

   bam.interval.list = list()
   for (bam.file in bam.files.c) {
      sample.name.str = file_path_sans_ext(basename(bam.file))
      for (interval in intervals.c) {
         temp = sprintf('%s(%s)', sample.name.str, str_extract(interval, pattern = '(?<=temp_)\\d+'))
         bam.interval.list[[temp]] = c(bam.file = bam.file, interval = interval)
      }
   }
   # bam.interval.list
   # $`123456-RN(0001)`
   # c(bam.file = 'bamfiles/1234567-RN.bam', interval = 'temp/germline.snp.RNA/temp_NNNN_of_N/scattered.interval_list')

   r = parallel::mclapply(bam.interval.list, FUN = ..RNAseqGermline.HaplotypeCaller, gatk = gatk, ref.fasta = ref.fasta, dbsnp = dbsnp, mc.cores = round(get_cpu_count()/3))
   message()
   print(r)
   message(r)
   if (!all(r == 0)) {
      quit(save = 'no', status = 1)
   }

   # output: {temp.folder}/{sample.name.str}/RNA.{NNNN}.filtered.germline.vcf.gz


   # ==========================Merge VCFs                              output: {output.folder}/{sample.name.str}.germline.vcf.gz========
   ..mergeVCFs <- function(sample.name.str) {

      vcf.files.c = list.files(path = sprintf('%s/%s', temp.folder, sample.name.str), pattern = '.filtered.germline.vcf.gz$', full.names = T)
      input.c = paste0(vcf.files.c, collapse = ' -I ')

      # merge========
      command.str = sprintf('%s --java-options -Xms4g MergeVcfs -I %s -O %s/%s.germline.vcf.gz', gatk, input.c, temp.folder, sample.name.str)
      message('\nRNA Germline MergeVCFs: ', sample.name.str)
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }

      # sort=========
      command.str = sprintf('%s --java-options -Xms4g SortVcf -I %s/%s.germline.vcf.gz -O %s/%s.germline.vcf.gz', gatk, temp.folder, sample.name.str, output.folder, sample.name.str)
      message('\nRNA Germline SortVcf: ', sample.name.str)
      r = system(command.str, ignore.stderr = T, ignore.stdout = T)
      if (r == 0) {
         message(sprintf('\n[%s] SUCCESS: %s', Sys.time(), command.str))
         return(0)
      } else {
         message(sprintf('\n[%s] Failed (%s): %s', Sys.time(), r, command.str))
         return(1)
      }
   }

   sample.name.list = sapply(file_path_sans_ext(basename(bam.files.c)), FUN = list)
   r = parallel::mclapply(sample.name.list, FUN = ..mergeVCFs, mc.cores = sample.each.time)
   message()
   print(r)
   message(r)
   if (!all(r == 0)) {
      quit(save = 'no', status = 1)
   }

}



#GATK.germline.CNV <- function() {}


# named.vector,
# strelka.path
# bed.file
# ref.fasta
# exome, manta.extra.arg = '', strelka.extra.arg = ''
Strelka.Somatic.Variants <- function(bams.list, parameter.file) {


   parameters.df = read.csv(parameter.file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)

   strelka.path = parameters.df['strelka.path', 'V2']
   bcftools.path = parameters.df['bcftools.path', 'V2']
   vcf2avinput.path = parameters.df['vcf2avinput.path', 'V2']

   ref.fasta = parameters.df['ref.fasta', 'V2']
   bed.file = parameters.df['bed.file', 'V2']

   is.exome = as.logical(parameters.df['is.exome', 'V2'])


   manta.extra.arg = parameters.df['manta.extra.arg', 'V2']
   strelka.extra.arg = parameters.df['strelka.extra.arg', 'V2']


   temp.folder = 'temp/somatic.snp.strelka'
   dir.create(temp.folder, recursive = T, showWarnings = F)

   output.folder = './somatic.snp.strelka'
   dir.create(output.folder, recursive = T, showWarnings = F)

   # read bams.list
   sample.list = generate.bam.list(bams.list)

   if (tools::file_ext(bed.file) == 'gz' & file.exists(sprintf('%s.tbi', bed.file))) {
      interval = bed.file
   } else {
      message('bed file must bgzipped and index by tabix.')
   }


   if (strelka.path == '') {
      configManta = 'configManta.py'
      configStrelka = 'configureStrelkaSomaticWorkflow.py'
   } else {
      configManta = sprintf('%s/configManta.py', strelka.path)
      configStrelka = sprintf('%s/configureStrelkaSomaticWorkflow.py', strelka.path)
   }

   r1 = system(sprintf('%s --version', configManta), ignore.stderr = T, ignore.stdout = T)
   r2 = system(sprintf('%s --version', configStrelka), ignore.stderr = T, ignore.stdout = T)
   if (r1 != 0 & r2 != 0) {
      message('Can not find Strelka')
      quit(save = 'no', status = 1)
   }

   # ..manta.and.strelka.somatic.workflows <- function(named.vector, strelka.path, bcftools.path, vcf2avinput.path, bed.file, ref.fasta, exome, manta.extra.arg = '', strelka.extra.arg = '')
   # strelka will use all cpus, parallel mode will crash it.
   for (name in names(sample.list)) {
      r = ..manta.and.strelka.somatic.workflows(named.vector = sample.list[[name]], strelka.path = strelka.path, bcftools.path = bcftools.path, vcf2avinput.path = vcf2avinput.path, bed.file = interval, ref.fasta = ref.fasta, exome = is.exome, manta.extra.arg = manta.extra.arg, strelka.extra.arg = strelka.extra.arg)
      message(sprintf('\n%s: %s', name, r))
   }

}





