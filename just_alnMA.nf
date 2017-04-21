import groovy.text.*
import java.io.*


params.aligners = "$baseDir/aligners.txt"
params.scores="$baseDir/scores.txt"
params.output_dir = ("$baseDir/output")
params.datasets_directory="$baseDir/benchmark_datasets"
datasets_home= file(params.datasets_directory)
params.score="bengen/baliscore"
params.bucket="100 200 300"
params.msa_dir = "MSAs"
params.msa_method ="msaprobs" //clustalo, mafft

//Reads which aligners to use
boxes = file(params.aligners).readLines().findAll { it.size()>0 }


//Reads which score to use
//boxes_score = file(params.scores).readLines().findAll { it.size()>0 }

boxes_score=["$params.score"]
/*
 * Creates a channel emitting a triple for each file in the datase composed
 * by the following element:
 *
 * -Name of the dataset (e.g. Balibase,Oxfam..)
 * -Name of the file (e.g. B11001_RV11 )
 * -the file itself
 *
 */

Channel
  .fromPath("${params.datasets_directory}/${params.dataset}/*.fa")
  .map { tuple( it.parent.name, it.baseName, it ) }
  .into{ dataset_fasta }

Channel
  .fromPath("${params.datasets_directory}/${params.dataset}/*.fa.ref")
  .map { tuple( it.name.tokenize('.')[0], it.parent.name, it ) }
  .into{ ref_fasta }

Channel.from( "${params.bucket}".tokenize() )
       .into{ bucs }
       .println()


bucs.combine(dataset_fasta)        
//    .println()
    .into{ data_pairs } 

/*
 * Execute an alignment job for each input sequence in the dataset
 */

process aln {

  tag "${params.msa_method}/${buc}/${id}"
  publishDir "${params.msa_dir}/$buc/$id"

  input:
  each method from boxes
  set buc, dataset_name, id, file(fasta) from data_pairs


  output:
  set id, buc, dataset_name, method, file('aln.fa') into alignments 
  file("*mafftdnd") 

  """

  #clustalo -i $fasta --outfmt fasta --guidetree-out the.tree --threads $task.cpus --force
  #t_coffee -seq $fasta -outfile aln.fa  -dpa -dpa_tree the.tree  -dpa_method ${params.msa_method}_msa -dpa_nseq $buc

  t_coffee -dpa -dpa_method ${params.msa_method} -dpa_tree mafftdnd -seq $fasta -outfile aln.fa  -dpa_nseq $buc
  """
}



ref_fasta.combine(alignments, by: 0)
//        .println()
	       .into{ data_aln}

/* 
 * Evaluate the alignment score 
 */

process score {
    tag "${params.score}/${buc}/${id}"
    container "${params.score}"
    publishDir "${params.msa_dir}/${buc}/${id}"
    
    input: 
    set ( id, loc, file(ref), buc, tag, method, file(aln) ) from data_aln
    

    output: 
    file("score.out")
    
    """
     msaa -r ${ref} -a ${aln} > score_temp.out
     cat score_temp.out | awk '{ print "SOP="\$2}' >  score.out
     
     perl ~/bin/extract_refAln.pl ${aln} ${ref} tmp_extracted.aln                          
     t_coffee -other_pg aln_compare -al1 ${ref} -al2 tmp_extracted.aln -compare_mode tc | grep -v "seq" | grep -v '*' | awk '{print "TC= "\$4}' >> score.out
    """
}

