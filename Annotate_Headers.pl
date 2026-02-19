#!/usr/bin/perl
use strict;
use warnings;

# Verifica os argumentos
if (@ARGV != 3) {
    die "Uso: perl script.pl <primeiro_multifasta> <segundo_multifasta_com_anotacoes> <saida_multifasta>\n";
}

my ($input_fasta, $annotated_fasta, $output_fasta) = @ARGV;

# Carrega as anotações SP+TM ou TM
my %annotations;
open(my $ann_fh, "<", $annotated_fasta) or die "Erro ao abrir $annotated_fasta: $!";
while (<$ann_fh>) {
    chomp;
    if (/^[^|]+\|([^|]+)\|.*\|\s*(SP\+TM|TM)$/) {
        $annotations{$1} = $2;
    }
}
close $ann_fh;

# Processa o primeiro multifasta
open(my $in_fh, "<", $input_fasta) or die "Erro ao abrir $input_fasta: $!";
open(my $out_fh, ">", $output_fasta) or die "Erro ao criar $output_fasta: $!";

while (<$in_fh>) {
    if (/^>([^|]+\|([^|]+)\|.*)/) {
        my $full_header = $1;
        my $id = $2;
        if (exists $annotations{$id}) {
            print $out_fh ">$full_header $annotations{$id}\n";
        } else {
            print $out_fh ">$full_header\n";
        }
    } else {
        print $out_fh $_;
    }
}

close $in_fh;
close $out_fh;
