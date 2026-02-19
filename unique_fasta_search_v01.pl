#!/usr/bin/perl
use strict;

#Globals
my $multifasta1 = $ARGV[0]; 	#multifasta path name
my $multifasta2 =  $ARGV[1];
my $how = $ARGV[2];  		# Flavor : seq, title or TFW  #  do you want to search unique based on sequences or titles?

				#>seq1 tyrosine receptor
				
				#TFW = title fist word, that will be >seq1

#usage
# script.pl infile1 infile2 flavor



my %fasta1;
my %fasta2;
my @unique1;
my @unique2;
my $before;
my @lines;
my @parts;
my $value ="";
my $seq = "garbage";
my $n = "0";
my $m1 = "0";
my $m2 = "0";

###reading input 
if ($how eq "seq"){
##multifasta1
    open (READ,$multifasta1) or die ("Error: $!");
	while (<READ>){
	    if ($_ =~ /^\>/){
		++$n;
		$fasta1{"$seq\n"}=$value;
		$value = $_;
		$seq ="";
	    }else{
		chomp $_;
		$seq = $seq.$_;
	    }
	}
    	$fasta1{"$seq\n"}=$value;
    close(READ);

    delete $fasta1{"garbage"};

##multifasta2
    open (READ,$multifasta2) or die ("Error: $!");
	while (<READ>){
	    if ($_ =~ /^\>/){
		++$n;
		$fasta2{"$seq\n"}=$value;
		$value = $_;
		$seq ="";
	    }else{
		chomp $_;
		$seq = $seq.$_;
	    }
	}
    	$fasta2{"$seq\n"}=$value;
    close(READ);

    delete $fasta2{"garbage"};


    ##finding uniques
    @unique1 = grep ! exists $fasta2{$_}, keys %fasta1;
    @unique2 = grep ! exists $fasta1{$_}, keys %fasta2;
    
    
    ##Printing final fasta sequences
    my $output = $multifasta1.".unique";
    open (WRITE,">$output") or die ("Error: $!");
    foreach (@unique1){
	next if ($_ =~ /garbage/);
	++$m1;
	print WRITE $fasta1{$_};
	print WRITE $_;
    }
    close (WRITE);
    
    my $output = $multifasta2.".unique";
    open (WRITE,">$output") or die ("Error: $!");
    foreach (@unique2){
	next if ($_ =~ /garbage/);
	++$m2;
	print WRITE $fasta2{$_};
	print WRITE $_;
    }
    
    ###Printing number of UNIQUE SEQUENCES
    print  "USING SEQUENCE COMPARISON\n$m1 unique sequences in fist file\n$m2 unique sequences in second file\n\n";
}


##reading input
if ($how eq "title"){
    open (READ,$multifasta1) or die ("Error: $!");
	while (<READ>){
	    if ($_ =~ /^\>/){
		++$n;
		$fasta1{$_}=$_;
		$before = $_;
	    }else{
		$fasta1{$before} = $fasta1{$before}.$_;
	    }
	}
    close(READ);


    open (READ,$multifasta2) or die ("Error: $!");
	while (<READ>){
	    if ($_ =~ /^\>/){
		++$n;
		$fasta2{$_}=$_;
		$before = $_;
	    }else{
		$fasta2{$before} = $fasta2{$before}.$_;
	    }
	}
    close(READ);



    ##finding uniques
    @unique1 = grep ! exists $fasta2{$_}, keys %fasta1;
    @unique2 = grep ! exists $fasta1{$_}, keys %fasta2;
    
    
    ##Printing final fasta sequences
    my $output = $multifasta1.".unique";
    open (WRITE,">$output") or die ("Error: $!");
    foreach (@unique1){
	++$m1;
	print WRITE $fasta1{$_};
    }
    close (WRITE);
    
    my $output = $multifasta2.".unique";
    open (WRITE,">$output") or die ("Error: $!");
    foreach (@unique2){
	++$m2;
	print WRITE $fasta2{$_};
    }
    
    ###Printing number of UNIQUE SEQUENCES
    print  "USING TITLE COMPARISON\n$m1 unique sequences in fist file\n$m2 unique sequences in second file\n\n";

}


##reading input
if ($how eq "TFW"){
    open (READ,$multifasta1) or die ("Error: $!");
	while (<READ>){
	    if ($_ =~ /^\>/){
		++$n;
		@parts = split (" ", $_);
		$parts[0] = $parts[0]."\n";
		#print $parts[0];
		$fasta1{$parts[0]}=$_;
		$before = $parts[0];
	    }else{
		$fasta1{$before} = $fasta1{$before}.$_;
	    }
	}
    close(READ);


    open (READ,$multifasta2) or die ("Error: $!");
	while (<READ>){
	    if ($_ =~ /^\>/){
		++$n;
		@parts = split (" ", $_);
		$parts[0] = $parts[0]."\n";
		#print $parts[0];
		$fasta2{$parts[0]}=$_;
		$before = $parts[0];
	    }else{
		$fasta2{$before} = $fasta2{$before}.$_;
	    }
	}
    close(READ);



    ##finding uniques
    @unique1 = grep ! exists $fasta2{$_}, keys %fasta1;
    @unique2 = grep ! exists $fasta1{$_}, keys %fasta2;
    
    
    ##Printing final fasta sequences
    my $output = $multifasta1.".unique";
    open (WRITE,">$output") or die ("Error: $!");
    foreach (@unique1){
	++$m1;
	print WRITE $fasta1{$_};
    }
    close (WRITE);
    
    my $output = $multifasta2.".unique";
    open (WRITE,">$output") or die ("Error: $!");
    foreach (@unique2){
	++$m2;
	print WRITE $fasta2{$_};
    }
    
    ###Printing number of UNIQUE SEQUENCES
    print  "USING TITLE FIRST WORD COMPARISON\n$m1 unique sequences in fist file\n$m2 unique sequences in second file\n\n";

}


exit;