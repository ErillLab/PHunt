#!/usr/bin/perl 

#Program to create PSSM and convert it's formats. 

use Data::Dumper;
use strict;

my $input = "@ARGV";

# Extracting input file name from command line

my @matrix;
if ($input =~ /-i\s*([^\s]+)/){
    my $infile = $1;
    open(MATRIX, $infile) or die "Cannot open file\"$infile\"\n\n ";
    @matrix = <MATRIX>;
    close MATRIX;

    unless (@matrix){
	print "The input file $infile is empty";	
	exit;
    }
} else {
    print "Specify input file name!\n";
    exit;
}

my $matrix = "@matrix";

@matrix = split(/[\s\n]*\/\/\/\/\/[\s\n]*/, $matrix);

#Extractin output file name from command line
my $outfile;
    if ($input =~ /-o\s*([^\s]+)/) {
	$outfile = $1;
    }

my $i;

if($matrix[0] eq "" && $matrix[1] ne ""){
shift(@matrix);
}


#Recognize the format of matrix, send the string containing 
#the matrix to appropriate subroutine, and store the results 
my $transfac_output;
my $fasta_output;
my $patser_output;
my $format;

for($i=0; $i<scalar(@matrix); $i++){
$matrix = $matrix[$i];

if($matrix ne " "){

if ($matrix =~ /^\s*ID.*\n/ ) {
$format = "tr";
    if($i==0){
    $transfac_output .= transfac_to_hash($matrix)."\n/////\n";;
    }else{
    $transfac_output .= transfac_to_hash($matrix);
    }

} elsif ($matrix =~ /\n*\s*(>.*\n)/) {
$format = "fa";
    my $name = $1;
    if($i==0){
    $fasta_output .= fasta_to_hash($matrix)."\n/////\n";
    }else{
    $fasta_output .= fasta_to_hash($matrix);
    }

} elsif ($matrix =~ /^\s*[Aa]\s*[:|]/) {
$format = "pa";
    if($i==0){
    $patser_output .= patser_to_promscan($matrix)."\n/////\n";
    }else{
    $patser_output .= patser_to_promscan($matrix);
    }

} else {
    print "Incorrect matrix format!\n";  
exit;
}    
}
}

open(OUTPUT, ">$outfile") or die  "Cannot open file \"$outfile\"\n\n";

if ($format eq "tr") {
print OUTPUT $transfac_output;
} elsif ($format eq "fa") {                             
print OUTPUT $fasta_output;
} elsif ($format eq "pa") {
print OUTPUT $patser_output;
}

close (OUTPUT);


###########################################################################################
#######################################
# Subroutines for convert
#######################################
###########################################################################################

# Definition of subroutine transfac_to_hash
#Subroutine transfac_to_hash takes the string with the original matrix in transfac format,
#creates the hash of arrays with values extracted from the matrix, and 
#sends them to next subroutine according to the user's choice  
sub transfac_to_hash {

    my ($matrixth) = @_;
        #Saving Name
        my $name;
        if ($matrixth =~ /(^\s*ID.*\n)/) {
    	    $name = $1;
	}
	
        #Begining fasta
        $name =~ s/ID/>/ig;
        $matrixth =~ s/ID.*\n//;
        $matrixth =~ s/BF.*\n//;
	$matrixth =~ s/[A-Za-z]*//g;
        $matrixth =~ s/\n\s*[0-9]{2}/\n/g;
        $matrixth =~ s/^\s*\n//g;
	
        # Counting length of matrix
        my $matrixlength;
        my $length;
        my $rows;
        my @matrix;
        $matrixlength = $matrix;
        $matrixlength =~ s/\///g;
	$matrixlength =~ s/ID.*\n//i;
        $matrixlength =~ s/BF.*\n//i;
        $matrixlength =~ s/\n\s[0-9]{2}//g;
	$matrixlength =~ s/[a-z]+//ig;
        $matrixlength =~ s/[0-9]+/$&,/g;
        $matrixlength =~ s/\n//g;
        $matrixlength =~ s/\s//g;
        $matrixlength =~ s/,/\n/g;
        @matrix = split /\s/, $matrixlength;
        $length = scalar @matrix;
        $rows = $length/4;
	
        #Initializing the arrays
    	my @am = "";
        my @cm = "";
	my @gm = "";
	my @tm = "";
    
			
	#Creating 2-D matrix
	my $i;
	for ($i = 0; $i < $rows;) {
      
	    $matrixth =~ /(\d+)\s*(\d+)\s*(\d+)\s*(\d+)/;
    	    push (@am, $1);
	    push (@cm, $2);
	    push (@gm, $3);
	    push (@tm, $4);
	    ++$i;
	    $matrixth =~ s/(\d+)\s*(\d+)\s*(\d+)\s*(\d+)//;
	}
	
	#Removing first (empty) character of each array
        shift @am;
        shift @cm;
        shift @gm;
        shift @tm;
			    
	#Forming hash of arrays
	my %matrixx;
	$matrixx{'a'} = [ @am ];
	$matrixx{'c'} = [ @cm ];
	$matrixx{'g'} = [ @gm ];
	$matrixx{'t'} = [ @tm ];

            $matrix = hash_to_promscan(%matrixx);
	    return $matrix;
}


#Definition of subroutine fasta_to_hash
#Subroutine fasta_to_hash takes the string with the original matrix or sequences in FASTA format,
#creates the hash of arrays with values extracted from the matrix (or sequences), and
#sends them to next subroutine according to the user's choice
sub fasta_to_hash {
    
    my ($matrixfh) = @_; 
    my $recognize = $matrixfh;
    $recognize =~ s/\s*>.*\n//;    
    #Initialization of hash of arrays
    my @a = "";
    my @c = "";
    my @g = "";
    my @t = "";

    #Recognizing the sequences
    if ($recognize =~ /\n*\s*[acgt]/i) {
	my $sequences = $matrixfh;
	my @sequences = split //, $sequences;
	
	#Counting the number of characters in matrices
	$sequences =~ /(>.*?)>/s;
	my $represent = $1;
	$represent =~ s/.*\n//;
	$represent =~ s/^\s*//;
	$represent =~ s/\s*$//;
	my @represent;
	@represent = split //, $represent;
	my $number_of_char = scalar @represent;

	#Counting the number of matrices

	my $i = 0;
	foreach (@sequences) {
	    if ( />/s ) {
	        ++$i;
	    }
	}

	my $number_of_matrices = $i;

	my $u;
	for ($u = 0 ; $u < $number_of_char; ++$u ) {
	    unshift (@a, "0");
	    unshift (@c, "0");
	    unshift (@g, "0");
	    unshift (@t, "0");
	}

	#Creating the hash of arrays
	my $e;

	$sequences = $sequences.">";
	for ($e = 0; $e < $number_of_matrices; ++$e) {
	    $sequences =~ /(>.*?)>/s;
	    my $one = $1;
	    $one =~ s/.*\n//;
	    $one =~ s/^\s*//;
	    $one =~ s/\s*$//;
	    my @one;
	    @one = split //, $one;
	
	    my $o = "0";
	    foreach (@one)  {
    		if ( /a/i ) {
	            ++$a[$o];
	        } elsif ( /c/i ) {
	            ++$c[$o];
	        } elsif ( /g/i ) {
	            ++$g[$o];
	        } elsif ( /t/i ) {
	            ++$t[$o];
	        }
	        ++$o;
	    }
	    $sequences =~ s/>.*?(>)/$1/s;
	}
    }



    else {
	
	#Handling the matrix
	$matrixfh =~ s/\s*>.*\n//;

        #Matching base wit appropriate row, deleting previous row,
        # and  removing whitespaces from the begining and the end of each line
        my $a;
        my $c;
        my $g;
        my $t;
        
        #a
        $matrixfh =~ /(\s*\d.*\n)/;
        $a = $1;
        $a =~ s/\s*$//;
        $a =~ s/^\s*//;
        @a = split /\s+/, $a;
        
        #c
        $matrixfh =~ s/(\s*\d.*\n)//;
        $matrixfh =~ /(\s*\d.*\n)/;
        $c = $1;
        $c =~ s/\s*$//;
        $c =~ s/^\s*//;
        @c = split /\s+/, $c;
        
	#g
        $matrixfh =~ s/(\s*\d.*\n)//;
        $matrixfh =~ /(\s*\d.*\n)/;
        $g = $1;
        $g =~ s/\s*$//;
        $g =~ s/^\s*//;
        @g = split /\s+/, $g;
    
	#t
	$matrixfh =~ s/.*\n//;
        $matrixfh =~ s/\s*$//;
        $matrixfh =~ s/^\s*//;
        @t = split /\s+/, $matrixfh;
	
    
    }
    #Forming the hash of arrays
    my %matrixx;
    $matrixx{'a'} = [ @a ];
    $matrixx{'c'} = [ @c ];
    $matrixx{'g'} = [ @g ];
    $matrixx{'t'} = [ @t ];


    #if ($format eq "pr") {
        $matrix = hash_to_promscan(%matrixx);
        return $matrix;
    #}
  
}

#Definition of subroutine patser
#Subroutine patser takes the string with the original matrix in patser or promscan format,
#creates the hash of arrays with values extracted from the matrix, and
#sends them to next subroutine if the user's choice is  Sequences, Transfac or Weblogo.
#Otherwise it just changes "|",":" or ">" sign.
sub patser_to_promscan {
    my ($matrix) = @_;

    #if ($format eq "pr") {
	$matrix =~ s/\|/:/g;
    #}

    return $matrix;

     
}


sub hash_to_promscan {
    my (%matrixx) = @_;
    unshift (@{ $matrixx{a} }, "A :");
    unshift (@{ $matrixx{c} }, "C :");
    unshift (@{ $matrixx{g} }, "G :");
    unshift (@{ $matrixx{t} }, "T :");
    $matrix = "@{ $matrixx{a} }\n@{ $matrixx{c} }\n@{ $matrixx{g} }\n@{ $matrixx{t} }";  
    return $matrix;
}

    								   			    	