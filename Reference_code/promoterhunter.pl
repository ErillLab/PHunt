#!/usr/bin/perl  -w

#####################################################################################################################
#                                                                          					    #	
#    PromoterHunter - A tool for promoter search in prokaryotic genomes      					    #	
#                                                                          					    #		
#####################################################################################################################
#                                                                                                                   #
#    This program is licensed under the Creative Commons Attribution-ShareAlike 3.0                                 # 
#    Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/         #
#    or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.     #
#                                                                                                                   #
#####################################################################################################################


use strict;
use Statistics::Descriptive;

##########################################################
#        	   Input processing                      #
##########################################################

my $input = "@ARGV";
my @seq_complete;
my @matrix;
my $seq_file;

if ($input =~ /-i\s*([^\s]+)/) {
    $seq_file = $1;
    open(FH, $seq_file) or die "Cannot open file \"$seq_file\"\n\n ";
    @seq_complete = <FH>;
    close FH;
}


my $seq_complete = "@seq_complete";
my $matrix_file;

if ($input =~ /-m\s*([^\s]+)/) {
    $matrix_file = $1;
    open(FH, $matrix_file) or die "Cannot open file \"$matrix_file\"\n\n ";
    @matrix = <FH>;
    close FH;
}

my $spacer_max;
my $spacer_min;
my $gc_content;
my $energy;
my $strand;
my $limit;

if ($input =~ /-x\s*([^-s\s]+)/) {
    $spacer_max = $1;
}

if ($input =~ /-n\s*([^-u\s]+)/) {
    $spacer_min = $1;
}

if ($input =~ /-g\s*([^-d\s]+)/) {
    $gc_content = $1;
}

if ($input =~ /-e\s*([^-l\s]+)/) {
    $energy = $1;
}

if ($input =~ /-s\s*([^-e\s]+)/) {
    $strand = $1;
}

if ($input =~ /-l\s*([^\s]+)/) {
    $limit = $1;
}

##########################################################
##			Warnings			##
##########################################################


if($seq_complete eq ""){
    print "No sequence to analyse!\n";
    exit;
}

unless($seq_complete =~ m/>/ ){
    print "Incorrect sequence format! Input sequence(s) has to be in FASTA format.\n";
    exit;
}


if($gc_content ne ""){
    if($gc_content >= 100){
	print "Incorrectly defined GC content!\n";
	exit;
    }unless($gc_content>0){
	print "Incorrectly defined GC content!\n";
	exit;
    }
}

if($spacer_max=~m/[^\d]/ || $spacer_max<1 || $spacer_max eq ""){
    print "Incorrectly defined distance between promoter regions!\n";
    exit;
}

if($spacer_min=~m/[^\d]/ || $spacer_min<1 || $spacer_min eq ""){
    print "Incorrectly defined distance between promoter regions!\n";
    exit;
}

if($spacer_max-$spacer_min>10 || $spacer_max-$spacer_min<0){
    print "Incorrectly defined distance between promoter regions!<br/>
    The difference between the upper and the lower limit of the distance between promoter regions has to be more than 0 and less than 11 bps.\n";
    exit;
}

if($limit ne "" && ($limit=~m/[^\d]/ || $limit<1 || $limit>100)){
    print "Incorrectly defined number of hits to be shown!\n";
    exit;
}

if($strand=~m/[^BS]/){
    print "Define whether both strands or single (positive) strand should be analysed [B/S].\n";
    exit;
}

if($energy=~ m/[^TF]/){
    print "Define whether Gibbs free energy should be calculated [T/F].\n";
    exit;
}

my @fasta_comments;
@fasta_comments = $seq_complete =~ />.*/g;

my @all_sequences;
@all_sequences = split(/>.*/, $seq_complete);
shift(@all_sequences);

my $total_length;

my $i=0;

foreach (@all_sequences){

$_ =~ s/\n//g;
$_ =~ s/\s//g;

if(length($_)>10000){
    print "Input sequence \"".$fasta_comments[$i]."\" is too long! The limit is 10000 letters per sequence.\n";
    exit;
}


if($_ =~ m/[^A|C|G|T|U|R|Y|S|W|K|M|B|D|H|V|N]/i){
    print "Input sequence \"".$fasta_comments[$i]."\" contains invalid letters!\n";
    exit;
}

$total_length += length($_);

$i++;
}

if($total_length>1000000){
    print "Too many data to be analyzed! The limit is 1000000 letters.\n";
    exit;
}

if(scalar(@all_sequences)>1000){
    print "Too many sequences to be analyzed! The limit is 1000 sequences.\n";
    exit;
}

system ("/var/www/cgi-bin/phisite/pssm-convert-hunter.pl -i ".$matrix_file."  -o ".$matrix_file."_f -f pr");

open(FH, $matrix_file."_f") or die "Cannot open file\"".$matrix_file."f\"\n\n ";
@matrix = <FH>;
close FH;

my $matrix1 = "@matrix";
my $matrix2 = "@matrix";
$matrix1 =~ s/\n*\/\/\/\/\/[\s\n\w\W]*//g;
$matrix2 =~ s/[\s\n\w\W]*\/\/\/\/\/\n*//g;


##########################################################
#        Initialization of variables                     #
##########################################################

my $j;
my $n;
my $m;
my $d;
my $p;

my $key;

my @IC1;
my @IC2;

my $ICad;
my $ICcy;
my $ICgu;
my $ICty;

my $num_of_seq;
my $motif_length1;
my $ad;
my $cy;
my $gu;
my $ty;
my @ad;
my @cy;
my @gu;
my @ty;
my @alignment_matrix1;

my $num_of_seq2;
my $motif_length2;
my $ad2;
my $cy2;
my $gu2;
my $ty2;
my @ad2;
my @cy2;
my @gu2;
my @ty2;
my @alignment_matrix2;

my %hash_of_total_scores;

my @array_of_all_energy_arrays;

my $shift;
my $run = 0;
my $complete_output;
my $stat;




foreach (@all_sequences){

undef %hash_of_total_scores;

my $seq = $_;
$seq =~ s/\n//g;
$seq =~ s/\s//g;

if ($seq =~ /[a-z]*/) {
$seq = uc $seq;
}

##########################################################
#        Accept ambiguous nucleotide symbols             #
##########################################################

if($seq =~ m/[^A|C|G|T]/){
my @R = ('A','G');
my @Y = ('C','T');
my @S = ('G','C');
my @W = ('A','T');
my @K = ('G','T');
my @M = ('A','C');
my @B = ('C','G','T');
my @D = ('A','G','T');
my @H = ('A','C','T');
my @V = ('A','C','G');
my @N = ('A','C','G','T');

my $R = $R[rand@R];
my $Y = $Y[rand@Y];
my $S = $S[rand@S];
my $W = $W[rand@W];
my $K = $K[rand@K];
my $M = $M[rand@M];
my $B = $B[rand@B];
my $D = $D[rand@D];
my $H = $H[rand@H];
my $V = $V[rand@V];
my $N = $N[rand@N];

$seq =~ tr/[U,R,Y,S,W,K,M,B,D,H,V,N]/[T,$R,$Y,$S,$W,$K,$M,$B,$D,$H,$V,$N]/;
}

if($gc_content eq ""){
    $gc_content = sprintf("%.1f",($seq =~ tr/C// + $seq =~ tr/G//)/length($seq)*100);
}


my $seq_length = length($seq);
my $out_sequence = $seq;


##########################################################
#        Calculate PWMs from input matrices              #
##########################################################

my $freq_a;
my $freq_t;
my $freq_c ;
my $freq_g;

if($gc_content>0){
    $freq_a = (1-($gc_content/100))/2;
    $freq_c = ($gc_content/100)/2;
    $freq_g = ($gc_content/100)/2;
    $freq_t = (1-($gc_content/100))/2;
}else{
    $freq_a = ($seq =~ tr/A//)/$seq_length;
    $freq_c = ($seq =~ tr/C//)/$seq_length;
    $freq_g = ($seq =~ tr/G//)/$seq_length;
    $freq_t = ($seq =~ tr/T//)/$seq_length;
}

$freq_a = sprintf("%.2f", $freq_a);
$freq_c = sprintf("%.2f", $freq_c);
$freq_g = sprintf("%.2f", $freq_g);
$freq_t = sprintf("%.2f", $freq_t);

$matrix1 =~ s/\s/\,/g; 
$matrix1 =~ s/\,+/\,/g;
my $num_matrix = $matrix1;
$num_matrix =~ s/[^\d |\,]//g; 
$num_matrix =~ s/^\,+//g;
my @matrix_array = split(',', $num_matrix);

$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@matrix_array);
my $sum_elements = sprintf("%.2f", $stat->sum());

$ad = $matrix1; 
$ad =~ s/C\,*[:|].+//g;
$ad =~ s/\,*A\,*[:|]//g;
$ad =~ s/^\,//g;
@ad = split(',', $ad);

$cy = $matrix1;
$cy =~ s/G\,*[:|].+//g;
$cy =~ s/\,*.+C\,*[:|]//g;
$cy =~ s/^\,//g;
@cy = split(',', $cy);

$gu = $matrix1;
$gu =~ s/T\,*[:|].+//g;
$gu =~ s/\,*.+G\,*[:|]//g;
$gu =~ s/^\,//g;
@gu = split(',', $gu);

$ty = $matrix1;
$ty =~ s/\,*.+T\,*[:|]//g;
$ty =~ s/^\,//g;
@ty = split(',', $ty);

@alignment_matrix1 = (@ad, @cy, @gu, @ty);

$motif_length1 = scalar(@ad);

if ($motif_length1 != 0){
    $num_of_seq = $sum_elements/$motif_length1;
}

if ($matrix2 ne ""){

    $matrix2 =~ s/\s/\,/g;
    $matrix2 =~ s/\,+/\,/g;
    my $num_matrix2 = $matrix2;
    $num_matrix2 =~ s/[^\d |\,]//g;
    $num_matrix2 =~ s/^\,+//g;
    my @matrix_array2 = split(',', $num_matrix2);

    $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@matrix_array2);
    my $sum_elements2 = sprintf("%.2f", $stat->sum());
						      						       
    $ad2 = $matrix2; 
    $ad2 =~ s/C\,*[:|].+//g;
    $ad2 =~ s/\,*A\,*[:|]//g;
    $ad2 =~ s/^\,//g;
    @ad2 = split(',', $ad2);

    $cy2 = $matrix2;
    $cy2 =~ s/G\,*[:|].+//g;
    $cy2 =~ s/\,*.+C\,*[:|]//g;
    $cy2 =~ s/^\,//g;
    @cy2 = split(',', $cy2);

    $gu2 = $matrix2;
    $gu2 =~ s/T\,*[:|].+//g;
    $gu2 =~ s/\,*.+G\,*[:|]//g;
    $gu2 =~ s/^\,//g;
    @gu2 = split(',', $gu2);

    $ty2 = $matrix2;
    $ty2 =~ s/\,*.+T\,*[:|]//g;
    $ty2 =~ s/^\,//g;
    @ty2 = split(',', $ty2);

    @alignment_matrix2 = (@ad2, @cy2, @gu2, @ty2);
    $motif_length2 = scalar(@ad2);
    if ($motif_length2 != 0){
	$num_of_seq2 = $sum_elements2/$motif_length2;
    }
}



##########################################################
#        Calculate information content of PWMs           #
##########################################################

for ($j = 0; $j < $motif_length1; $j++) {

    if ($ad[$j] != 0 && $freq_a != 0){
	$ICad = ($ad[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j])) * log(($ad[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j]))/$freq_a);
    }else{
	$ICad = 0;
    }

    if ($cy[$j] != 0 && $freq_c != 0){
	$ICcy = ($cy[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j])) * log(($cy[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j]))/$freq_c);
    }else{
	$ICcy = 0;
    }

    if ($gu[$j] != 0 && $freq_g != 0){
	$ICgu = ($gu[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j])) * log(($gu[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j]))/$freq_g);
    }else{
	$ICgu = 0;
    }

    if ($ty[$j] != 0 && $freq_t != 0){
	$ICty =  ($ty[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j])) * log(($ty[$j]/($ad[$j]+$cy[$j]+$gu[$j]+$ty[$j]))/$freq_t);
    }else{
	$ICty = 0;
    }

    my $IC_sum = $ICad + $ICcy + $ICgu + $ICty;
    push (@IC1, $IC_sum );
}

if ($matrix2 ne ""){
    for ($j = 0; $j < $motif_length2; $j++) {
	if ($ad2[$j] != 0 && $freq_a != 0){
	    $ICad = ($ad2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j])) * log(($ad2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j]))/$freq_a);
	}else{
	    $ICad = 0;
	}	

	if ($cy2[$j] != 0 && $freq_c != 0){
	    $ICcy = ($cy2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j])) * log(($cy2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j]))/$freq_c);
	}else{
	    $ICcy = 0;
	}

	if ($gu2[$j] != 0 && $freq_g != 0){
	    $ICgu = ($gu2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j])) * log(($gu2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j]))/$freq_g);
	}else{
	    $ICgu = 0;
	}

	if ($ty2[$j] != 0 && $freq_t != 0){
	    $ICty =  ($ty2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j])) * log(($ty2[$j]/($ad2[$j]+$cy2[$j]+$gu2[$j]+$ty2[$j]))/$freq_t);
	}else{
	    $ICty = 0;
	}

	my $IC_sum = $ICad + $ICcy + $ICgu + $ICty;
	push (@IC2, $IC_sum );
    }
}

$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@IC1);
my $IC1 = sprintf("%.2f", $stat->sum());
my $IC2;

if ($matrix2 ne ""){
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@IC2);
$IC2 = sprintf("%.2f", $stat->sum());
}


##########################################################
#        Calculate weights in PWMses                     #
##########################################################

if ($freq_a != 0) {
 foreach (@ad) {
 $_ = (($_ + $freq_a) / ($num_of_seq+1))/$freq_a;
 $_ = log($_);
	$_ = sprintf("%.2f", $_);
 }
 if ($matrix2 ne ""){
 foreach (@ad2) {
 $_ = (($_ + $freq_a) / ($num_of_seq2+1))/$freq_a;
 $_ = log($_);
	$_ = sprintf("%.2f", $_);
	
 }}
}

if ($freq_c != 0) {
 foreach (@cy) {
  $_ = (($_ + $freq_c) / ($num_of_seq+1))/$freq_c;
   $_ = log($_); 
	$_ = sprintf("%.2f", $_);
 }
 if ($matrix2 ne ""){
 foreach (@cy2) {
  $_ = (($_ + $freq_c) / ($num_of_seq2+1))/$freq_c;
   $_ = log($_); 
	$_ = sprintf("%.2f", $_);
    }}
}

if ($freq_g != 0) {    
 foreach (@gu) {
  $_ = (($_ + $freq_g) / ($num_of_seq+1))/$freq_g;
   $_ = log($_);
	$_ = sprintf("%.2f", $_);
 }
 if ($matrix2 ne ""){
 foreach (@gu2) {
   $_ = (($_ + $freq_g) / ($num_of_seq2+1))/$freq_g;
    $_ = log($_);
	$_ = sprintf("%.2f", $_);
    }}
}
     
if ($freq_t != 0) {     
 foreach (@ty) {
  $_ = (($_ + $freq_t) / ($num_of_seq+1))/$freq_t;
   $_ = log($_);
	$_ = sprintf("%.2f", $_);
 }
 if ($matrix2 ne ""){
 foreach (@ty2) {
  $_ = (($_ + $freq_t) / ($num_of_seq2+1))/$freq_t;
   $_ = log($_);
	$_ = sprintf("%.2f", $_);
	
    }}
}


my @ad_rev = reverse(@ad);
my @cy_rev = reverse(@cy);
my @gu_rev = reverse(@gu);
my @ty_rev = reverse(@ty); 
my @ad2_rev = reverse(@ad2);
my @cy2_rev = reverse(@cy2);
my @gu2_rev = reverse(@gu2);
my @ty2_rev = reverse(@ty2);

my @max_array1;
my @min_array1;
my $max1 = 0;
my $min1 = 0;

my @max_array2;
my @min_array2;
my $max2 = 0;
my $min2 = 0;

my $mine=0;
my $maxe=0;

my $sig_level;

my $score_median;
my $score_sd;

my $threshold_1;
my $threshold_2;


##########################################################
#        Calculate minimal and maximal PWM scores        #
##########################################################

for ($i = 0; $i < $motif_length1; $i++) {
my @column_weights;
push (@column_weights, $ad[$i], $cy[$i], $gu[$i], $ty[$i]); 
@column_weights = sort{$b <=> $a}(@column_weights); 
push (@max_array1, $column_weights[0]);
@column_weights = sort{$a <=> $b}(@column_weights);
push (@min_array1, $column_weights[0]);
}
 
 foreach (@max_array1) {
   $max1 += $_;
   }
   
   foreach (@min_array1) {
      $min1 += $_;
         }


if($matrix2 ne ""){
for ($i = 0; $i < $motif_length2; $i++) {
my @column_weights2;
push (@column_weights2, $ad2[$i], $cy2[$i], $gu2[$i], $ty2[$i]);
@column_weights2 = sort{$b <=> $a}(@column_weights2);
push (@max_array2, $column_weights2[0]);
@column_weights2 = sort{$a <=> $b}(@column_weights2);
push (@min_array2, $column_weights2[0]);
}

 foreach (@max_array2) {
    $max2 += $_;
    }
       
    foreach (@min_array2) {
	$min2 += $_;
	    }
}


my @seq = split ('', $seq);


##########################################################
#        Initialization of variables                     #
##########################################################

my @array_of_subseqs;
my @array_of_subseqs2;

my @array_of_weights;
my @array_of_weights_rev;
my @array_of_weights2;
my @array_of_weights2_rev;

my @array_of_sums_of_weights;
my @array_of_sums_of_weights_rev;
my @array_of_sums_of_weights2;
my @array_of_sums_of_weights2_rev;
my @begin_array;
my @begin_array_rev;
my @begin_array2;
my @begin_array2_rev;
my @strand;
my @strand_rev;
my @strand2;
my @strand2_rev;

my @cutoff_begins;
my @cutoff_begins2;
my @cutoff_begins_rev;
my @cutoff_begins2_rev;
my @cutoff_scores;
my @cutoff_scores2;
my @cutoff_scores_rev;
my @cutoff_scores2_rev;
my @cutoff_strand;
my @cutoff_strand_rev;

my @summary_score;
my @summary_score_rev;
my @score1;
my @score2;

my @norm_summary_score;
my @norm_energy;

my @final_score;

my @energy_array;
my @resulting_energy;

my @all_weights;
my @all_weights2;


my @scScoresM;
my @scScoresM2;
my @scScoresM_rev;
my @scScoresM2_rev;
my @scScores1;
my @scScores2;



my $w = 0;

##########################################################
#            Calculate Gibbs free energy                 #
##########################################################
if($energy eq "T"){
$w = 50;

#Double each character:
    $seq =~ s/A/AA/g;
    $seq =~ s/C/CC/g;
    $seq =~ s/G/GG/g;
    $seq =~ s/T/TT/g;

#Delete the first character of the string:
    $seq =~ s/^\w//;

#Delete the last character of the string:
    $seq =~ s/\w$//;

#Create two-character groups:
    $seq =~ s/\w{2}/$& /g;

#Replace with appropriate free-energy values:
    $seq =~ s/AA/4.18/g;
    $seq =~ s/AC/6.02/g;
    $seq =~ s/AG/5.36/g;
    $seq =~ s/AT/3.68/g;
    $seq =~ s/CA/6.07/g;
    $seq =~ s/CC/7.7/g;
    $seq =~ s/CG/9.08/g;
    $seq =~ s/CT/5.36/g;
    $seq =~ s/GA/5.44/g;
    $seq =~ s/GC/9.37/g;
    $seq =~ s/GG/7.7/g;
    $seq =~ s/GT/6.02/g;
    $seq =~ s/TA/2.43/g;
    $seq =~ s/TC/5.44/g;
    $seq =~ s/TG/6.07/g;
    $seq =~ s/TT/4.18/g;

#Convert string to array:
    my @seq_array;
    @seq_array = split(' ', $seq);

#Counting free energy within sliding-window:

    for($i=0; $i<($seq_length-$w); $i++) {

	my @seq_win;
        @seq_win = @seq_array[$i..($w+$i-1)];
	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@seq_win);
	my $energy = sprintf("%.2f", $stat->mean());

        push (@energy_array, $energy);
    }

    my $position;
    if($w%2==0){
	$position = $w/2;
    }else{
	$position = ($w+1)/2;
    }

    $array_of_all_energy_arrays[$run] = [@energy_array];

    $maxe = 3;
    $mine = 9.2;

}

##########################################################
#        Scan sequence with the first PWM	         #
##########################################################



for ($i = 0; $i < (scalar @seq - $motif_length1 + 1); $i++) {

    my @subseq;
    @subseq= @seq[$i..($motif_length1 + $i - 1)];    

    push (@array_of_subseqs, "@subseq");

	for ($m = 0; $m < $motif_length1; $m++) {
	    if ($strand eq "B"){    
		if ($subseq[$m] eq "A") {
	    
		    push (@array_of_weights, $ad[$m]);
		    push (@array_of_weights_rev, $ty_rev[$m]);
	    
		}elsif ($subseq[$m] eq "C") {
	    
		    push (@array_of_weights, $cy[$m]);
		    push (@array_of_weights_rev, $gu_rev[$m]);
	    
		}elsif ($subseq[$m] eq "G") {
	    
		    push (@array_of_weights, $gu[$m]);
		    push (@array_of_weights_rev, $cy_rev[$m]);
	    
		}else{ 

		    push (@array_of_weights, $ty[$m]);
		    push (@array_of_weights_rev, $ad_rev[$m]);
	    
		}
	    }else{

		if ($subseq[$m] eq "A") {

		    push (@array_of_weights, $ad[$m]);

		}elsif ($subseq[$m] eq "C") {

		    push (@array_of_weights, $cy[$m]);

		}elsif ($subseq[$m] eq "G") {

		    push (@array_of_weights, $gu[$m]);

		}else{

		    push (@array_of_weights, $ty[$m]);

		}
	    }	     		
	}

    my $sum_of_weights = 0;
    foreach (@array_of_weights) {
	$sum_of_weights += $_;
        $sum_of_weights = sprintf("%.2f", $sum_of_weights);
    }
   
    my $sum_of_weights_rev = 0;
    if($strand eq "B"){
        foreach (@array_of_weights_rev) {
    	    $sum_of_weights_rev += $_;
	    $sum_of_weights_rev = sprintf("%.2f", $sum_of_weights_rev);     
	}    
    }
					   
    @array_of_weights = 0;

    if($strand eq "B"){
	@array_of_weights_rev = 0;
    }


    push(@all_weights, $sum_of_weights);
    if($strand eq "B"){
	push(@all_weights, $sum_of_weights_rev);
    }

    push (@array_of_sums_of_weights, $sum_of_weights);

    if($strand eq "B"){
	push (@array_of_sums_of_weights_rev, $sum_of_weights_rev);
    }

}

$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@all_weights);
$score_sd = sprintf("%.2f", $stat->standard_deviation());
$score_median = sprintf("%.2f", $stat->median());
$threshold_1 = $score_median + $score_sd;

my @tmp_array;
my @tmp_array_rev;

for($i=0; $i<scalar(@array_of_sums_of_weights); $i++){
    if($array_of_sums_of_weights[$i]>$threshold_1){
	push (@begin_array, $i+1);
	push (@strand, "+");
	push (@tmp_array,$array_of_sums_of_weights[$i] );
    }
}

for($i=0; $i<scalar(@array_of_sums_of_weights_rev); $i++){
    if($array_of_sums_of_weights_rev[$i]>$threshold_1){
	push (@begin_array_rev, $i+1);
	push (@strand_rev, "-");
	push (@tmp_array_rev,$array_of_sums_of_weights_rev[$i] );
    }
}

@array_of_sums_of_weights = @tmp_array;
@array_of_sums_of_weights_rev = @tmp_array_rev;

if($matrix2 eq ""){

    for ($i=0; $i < scalar @array_of_sums_of_weights; $i++){
	push(@scScoresM, sprintf("%.2f",($array_of_sums_of_weights[$i]-$min1)/($max1-$min1)));
    }

    if($strand eq "B"){
	for ($i=0; $i < scalar @array_of_sums_of_weights_rev; $i++){
	    push(@scScoresM_rev, sprintf("%.2f",($array_of_sums_of_weights_rev[$i]-$min1)/($max1-$min1)));
	}
    }


    if($strand eq "B"){
	@summary_score = (@array_of_sums_of_weights, @array_of_sums_of_weights_rev);
        @begin_array2 = (@begin_array, @begin_array_rev);
	@strand = (@strand, @strand_rev);
	@scScores1 = (@scScoresM, @scScoresM_rev);
    }else{
	@summary_score = @array_of_sums_of_weights;
        @begin_array2 = @begin_array;
	@strand = @strand;
        @scScores1 = @scScoresM;
    }
}


##########################################################
#        Scan sequence with the second PWM               #
##########################################################

if ($matrix2 ne ""){
for ($i = 0; $i < (scalar @seq - $motif_length2 + 1); $i++) {

    my @subseq;
	@subseq= @seq[$i..($motif_length2 + $i - 1)];
	 
	push (@array_of_subseqs2, "@subseq");
	 
	for ($m = 0; $m < $motif_length2; $m++) {
	    if ($strand eq "B"){
		if ($subseq[$m] eq "A") {
		    push (@array_of_weights2, $ad2[$m]);
		    push (@array_of_weights2_rev, $ty2_rev[$m]);
	 
		}elsif ($subseq[$m] eq "C") {
		    push (@array_of_weights2, $cy2[$m]);
		    push (@array_of_weights2_rev, $gu2_rev[$m]);
	
		}elsif ($subseq[$m] eq "G") {
		    push (@array_of_weights2, $gu2[$m]);
		    push (@array_of_weights2_rev, $cy2_rev[$m]);
	 
		}else{ 
		    push (@array_of_weights2, $ty2[$m]); 
		    push (@array_of_weights2_rev, $ad2_rev[$m]); 

		}
	    }else{
         
		if ($subseq[$m] eq "A") {
		    push (@array_of_weights2, $ad2[$m]);

		}elsif ($subseq[$m] eq "C") {
		    push (@array_of_weights2, $cy2[$m]);

		}elsif ($subseq[$m] eq "G") {
		    push (@array_of_weights2, $gu2[$m]);

		}else{
        	    push (@array_of_weights2, $ty2[$m]);
		}
     	    }	 
	}


    my $sum_of_weights = 0;
    foreach (@array_of_weights2) {
	$sum_of_weights += $_;
        $sum_of_weights = sprintf("%.2f", $sum_of_weights);
    }

    my $sum_of_weights_rev = 0;
    if($strand eq "B"){
	foreach (@array_of_weights2_rev) {
	    $sum_of_weights_rev += $_;
	    $sum_of_weights_rev = sprintf("%.2f", $sum_of_weights_rev);
	}
    }										     

    @array_of_weights2 = 0;

    if($strand eq "B"){
	@array_of_weights2_rev = 0;
    }

    push(@all_weights2, $sum_of_weights);

    if($strand eq "B"){
	push(@all_weights2, $sum_of_weights_rev);
    }
  
    push (@array_of_sums_of_weights2, $sum_of_weights);

    if($strand eq "B"){
	push (@array_of_sums_of_weights2_rev, $sum_of_weights_rev);
    }
}										     

    $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@all_weights2);
    $score_sd = sprintf("%.2f", $stat->standard_deviation());
    $score_median = sprintf("%.2f", $stat->median());
    $threshold_2 = $score_median+$score_sd;

    my @tmp_array2;
    my @tmp_array2_rev;

    for($i=0; $i<scalar(@array_of_sums_of_weights2); $i++){
	if($array_of_sums_of_weights2[$i]>$threshold_2){
	    push (@begin_array2, $i+1);
	    push (@strand2, "+");
	    push (@tmp_array2,$array_of_sums_of_weights2[$i] );
	}
    }

    for($i=0; $i<scalar(@array_of_sums_of_weights2_rev); $i++){
	if($array_of_sums_of_weights2_rev[$i]>$threshold_2){
	    push (@begin_array2_rev, $i+1);
	    push (@strand2_rev, "-");
	    push (@tmp_array2_rev,$array_of_sums_of_weights2_rev[$i] );
	}
    }

    @array_of_sums_of_weights2 = @tmp_array2;
    @array_of_sums_of_weights2_rev = @tmp_array2_rev;


#############################################################################################
#    Choose only those sites(begins and scores), which have appropriate mutual distance     #
#############################################################################################

    my $temp;
    for ($i=0; $i < scalar @array_of_sums_of_weights; $i++){
	for ($n=$spacer_min; $n<=$spacer_max; $n++) {
	    if (grep {$_ eq $begin_array[$i]+$motif_length1+$n} @begin_array2) {
		push (@cutoff_begins, $begin_array[$i]);
		push (@cutoff_begins2, $begin_array[$i]+$motif_length1+$n);
		my $elm;
		    for ($temp=0; $temp <scalar @begin_array2; $temp++){
			if ($begin_array[$i]+$motif_length1+$n == $begin_array2[$temp]){
			    $elm = $temp;
			}
		    }

		push (@cutoff_scores, $array_of_sums_of_weights[$i]);
		push (@cutoff_scores2, $array_of_sums_of_weights2[$elm]);
		push (@scScoresM, sprintf("%.2f",($array_of_sums_of_weights[$i]-$min1)/($max1-$min1)));
		push (@scScoresM2, sprintf("%.2f",($array_of_sums_of_weights2[$elm]-$min2)/($max2-$min2)));
		push (@cutoff_strand, "+");
	    }
	}
    }


    if($strand eq "B"){
	for ($i=0; $i < scalar @array_of_sums_of_weights2_rev; $i++){
	    for ($n=$spacer_min; $n<=$spacer_max; $n++) {

		if (grep {$_ eq $begin_array2_rev[$i]+$motif_length2+$n} @begin_array_rev) {
		    push (@cutoff_begins_rev , $begin_array2_rev[$i]+$n+$motif_length2);
		    push (@cutoff_begins2_rev, $begin_array2_rev[$i]);
		    my $elm;
			for ($temp=0; $temp <scalar @begin_array_rev; $temp++){
			    if ($begin_array2_rev[$i]+$motif_length2+$n == $begin_array_rev[$temp]){
				$elm = $temp;
			    }
			}

		    push (@cutoff_scores_rev, $array_of_sums_of_weights_rev[$elm]);
		    push (@cutoff_scores2_rev, $array_of_sums_of_weights2_rev[$i]);
		    push (@scScoresM_rev, sprintf("%.2f",($array_of_sums_of_weights_rev[$elm]-$min1)/($max1-$min1)));
		    push (@scScoresM2_rev, sprintf("%.2f",($array_of_sums_of_weights2_rev[$i]-$min2)/($max2-$min2)));
		    push (@cutoff_strand_rev, "-");
		}
	    }
	}
    }




#############################################################
#            Sum scores (PWM1+PWM2)			    #	
#############################################################


    for($i=0; $i < scalar@cutoff_scores; $i++){
	my $summary_score = $cutoff_scores[$i]+$cutoff_scores2[$i];
	$summary_score = sprintf("%.2f", $summary_score);
	push(@summary_score, $summary_score);
    }

    if($strand eq "B"){
	for($i=0; $i < scalar@cutoff_scores2_rev; $i++){
	    my $summary_score = $cutoff_scores2_rev[$i]+$cutoff_scores_rev[$i];
	    $summary_score = sprintf("%.2f", $summary_score);
	    push(@summary_score_rev, $summary_score);
	}
    }
#############################################################
#            Join direct and reverse arrays                 #
#############################################################


    if($strand eq "B"){
    
	@summary_score = (@summary_score, @summary_score_rev);
        @begin_array = (@cutoff_begins, @cutoff_begins_rev);
        @begin_array2 = (@cutoff_begins2, @cutoff_begins2_rev);
	@score1 = (@cutoff_scores, @cutoff_scores_rev);
	@score2 = (@cutoff_scores2, @cutoff_scores2_rev);
	@strand = (@cutoff_strand, @cutoff_strand_rev);
        @scScores1 = (@scScoresM,@scScoresM_rev);
        @scScores2 = (@scScoresM2,@scScoresM2_rev);

    }else{
	
	@summary_score = @summary_score;
	@begin_array = @cutoff_begins;
        @begin_array2 = @cutoff_begins2;
	@score1 = @cutoff_scores;
	@score2 = @cutoff_scores2;
        @strand = @cutoff_strand;
        @scScores1 = @scScoresM;
        @scScores2 = @scScoresM2;

    }

}		  



#############################################################
# Calculate mean free energy round the potential promoters  #
#############################################################

if($w%2==0){
    $shift = $w/2;
}else{
    $shift = ($w+1)/2;
}


my $k = 0;
    foreach (@begin_array2){
	
	if($energy eq "T"){

	    my $mean_energy;
	    if($strand[$k] eq "+"){

		my $left = $_-100-$shift+$motif_length2;
		if($left<0){$left = 0;} 
		my $right = $_+30-$shift+$motif_length2;
		if($right>scalar(@energy_array)-1){$right = scalar(@energy_array)-1;}
		my @energy_region = @energy_array[$left..$right];
	        $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@energy_region);
            	$mean_energy = sprintf("%.2f", $stat->mean());

	    }else{

		my $left = $_-30-$shift;
		if($left<0){$left = 0;}
		my $right = $_+100-$shift;
		if($right>scalar(@energy_array)-1){$right = scalar(@energy_array)-1;}
		my @energy_region = @energy_array[$left..$right];
                $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@energy_region);
                $mean_energy = sprintf("%.2f", $stat->mean());

	    }

        push(@resulting_energy, $mean_energy);
	$k++;

    }else{
	push(@resulting_energy, 0);

    }

}

#############################################################
#               Calculate final scores			    #
#############################################################

if($matrix2 ne "" && scalar(@scScores1)>0 ){
    for($i=0; $i<scalar(@scScores1); $i++){
	push (@norm_summary_score, ($scScores1[$i]+$scScores2[$i]));
    }
}else{
    @norm_summary_score = @scScores1;
}


#############################################################
#        Calculate final scores (PWM + free energy)	    #
#############################################################

if($energy eq "T"){

    if(scalar(@resulting_energy)>0){
        foreach(@resulting_energy){
	    if ($_ ne ""){
	        my $norm = sprintf("%.3f",($_ - $mine)/($maxe-$mine));
	        push (@norm_energy, $norm);

	    }
	}
    }

    for($i=0; $i<scalar @norm_summary_score; $i++){
	push(@final_score, sprintf("%.2f",($norm_summary_score[$i]+2*$norm_energy[$i] )));
    }


}else{
	for($i=0; $i<scalar @norm_summary_score; $i++){
	    push(@final_score, sprintf("%.2f",($norm_summary_score[$i])));
	}
}


#############################################################
#        	Set significance level                      #
#############################################################

my $el = 0.6;
my $ml = 0.8;

if($matrix2 eq ""){
    if($energy eq "T"){
	$sig_level = $ml+(($el+0.5-$gc_content/100)*2);
    }else{
	$sig_level = $ml;
    }
}else{
    if($energy eq "T"){
        $sig_level = 2*$ml+(($el+0.5-$gc_content/100)*2); 
    }else{
	$sig_level = 2*$ml;
    }
}

#############################################################
#             Put results to hash and order them            #
#############################################################
#			single PWM scan                     #
#############################################################


    my $csv;
    my $key;
    my $l;
    my $sig_mark;

if ($matrix2 eq "") {
    for($d=0; $d < scalar @summary_score; $d++){
	$hash_of_total_scores{$d} = ({'begin' => $begin_array2[$d], 'strand' => $strand[$d], 'summary_score' => $summary_score[$d], 'energy'=>$resulting_energy[$d], 'final_score' => $final_score[$d]});
    }

    $csv = "";
    $l=0;

    foreach $key (sort order keys %hash_of_total_scores) {

	if($hash_of_total_scores{$key}->{'final_score'}>=$sig_level){
	    $sig_mark = "!";
	}else{
	    $sig_mark = "";
	}

	if($l<$limit || $limit eq ""){
	    $csv .= $hash_of_total_scores{$key}->{'strand'}.",".$hash_of_total_scores{$key}->{'begin'}.",".$hash_of_total_scores{$key}->{'summary_score'}.",".$hash_of_total_scores{$key}->{'energy'}.",".$hash_of_total_scores{$key}->{'final_score'}.",".$sig_mark.",\n";
	    $l++;
	}
    }

    $complete_output .= "#ID:$fasta_comments[$run]\n#Sequence:$out_sequence\n#Seq_length:$seq_length\n#Matrix1:@alignment_matrix1\n#PWM1:@ad @cy @gu @ty\n#Min1:$min1\n#Max1:$max1\n#IC1:$IC1\n#Motif_length1:$motif_length1\nStrand,Begin,Score,Energy,Final score\n$csv";


#############################################################
#                       double PWM scan                     #
#############################################################

}else{

    for ($d=0; $d < scalar @summary_score; $d++) {
	$hash_of_total_scores{$d} = ({ 'begin1' => $begin_array[$d], 'begin2' => $begin_array2[$d], 'strand' => $strand[$d], 'score1' => $score1[$d], 'score2' => $score2[$d], 'summary_score' => $summary_score[$d], 'energy' => $resulting_energy[$d], 'final_score' => $final_score[$d]}); 
    }

    $csv = "";
    $l = 0;
    my $push;
    my @existingBegins1;
    my @existingBegins2;

    foreach $key (sort order keys %hash_of_total_scores) {
	$push = 1;

    if($hash_of_total_scores{$key}->{'begin2'}){
	if(grep $_ eq $hash_of_total_scores{$key}->{'begin2'}, @existingBegins2){
	    $push = 0;
	}
    }
    if($hash_of_total_scores{$key}->{'begin1'}){

        if(grep $_ eq $hash_of_total_scores{$key}->{'begin1'}, @existingBegins1 ){
	    $push = 0;
	}
    }
        push(@existingBegins2, $hash_of_total_scores{$key}->{'begin2'});
	push(@existingBegins1, $hash_of_total_scores{$key}->{'begin1'});

	if($hash_of_total_scores{$key}->{'final_score'}>=$sig_level){
	    $sig_mark = "!";
	}else{
	    $sig_mark = "";
	}

	if(($push==1 && $l<$limit) || ($push==1 && $limit eq "")){
	    $csv .= $hash_of_total_scores{$key}->{'strand'}.",".$hash_of_total_scores{$key}->{'begin1'}.",".$hash_of_total_scores{$key}->{'score1'}.",".$hash_of_total_scores{$key}->{'begin2'}.",".$hash_of_total_scores{$key}->{'score2'}.",".$hash_of_total_scores{$key}->{'summary_score'}.",".$hash_of_total_scores{$key}->{'energy'}.",".$hash_of_total_scores{$key}->{'final_score'}.",".$sig_mark.",\n";
	    $l++;

	}
    }

    $complete_output .= "#ID:$fasta_comments[$run]\n#Sequence:$out_sequence\n#Seq_length:$seq_length\n#Matrix1:@alignment_matrix1\n#Matrix2:@alignment_matrix2\n#PWM1:@ad @cy @gu @ty\n#PWM2:@ad2 @cy2 @gu2 @ty2\n#Min1:$min1\n#Max1:$max1\n#Min2:$min2\n#Max2:$max2\n#IC1:$IC1\n#IC2:$IC2\n#Motif_length1:$motif_length1\n#Motif_length2:$motif_length2\nStrand,Begin 1,Score 1,Begin 2,Score 2,Summary score,Energy score,Final score\n$csv";
}


$run++;


}


sub order {
  $hash_of_total_scores{$b}->{final_score} <=> $hash_of_total_scores{$a}->{final_score};
}

####################################################
#       		Output                     #
####################################################

if ($input =~ /-o\s*([^-g\s]+)/) {
    my  $outfile = $1;
    open(FH, ">$outfile") or die "can't open $outfile: $!\n";
    print FH $complete_output;
    close (FH);
} else {
    die "Specify output file name";
}


if($energy eq "T"){

    my @sorted_array_of_all_energy_arrays = sort{ scalar @{$b} <=>  scalar @{$a}} (@array_of_all_energy_arrays);
    my $energy_csv = "";

    foreach(@fasta_comments){
	$energy_csv .= $_.":position,".$_.":freeEnergy,";
    }

    $energy_csv .= "\n";

    for($i=0; $i<scalar(@{$sorted_array_of_all_energy_arrays[0]}); $i++){
	for($p=0; $p<scalar(@fasta_comments); $p++){
	    if(exists($array_of_all_energy_arrays[$p][$i])){
		$energy_csv .= ($shift+$i).",".$array_of_all_energy_arrays[$p][$i].",";
	    }else{
		$energy_csv .= ",,";
	    }
	}
	$energy_csv .= "\n";
    }

    if ($input =~ /-o\s*([^-g\s]+)/) {

	my  $outfile_energy = $1."_e";
        open(FH, ">$outfile_energy") or die "can't open $outfile_energy: $!\n";
	print FH $energy_csv;
        close (FH);
    } else {
	die "Specify output file name";
    }

}


