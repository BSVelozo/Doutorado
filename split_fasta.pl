#!/usr/bin/perl
use strict;
use warnings;

my $input_file = 'gencode.fasta';  # Substitua pelo nome do seu arquivo multifasta

open my $fh, '<', $input_file or die "Não foi possível abrir o arquivo '$input_file': $!";

my $output_file;
while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^>(.*)/) {  # Verifica se a linha começa com '>'
        if ($output_file) {
            close OUT;  # Fecha o arquivo anterior se houver
        }
        $output_file = $1 . '.fasta';  # Define o novo nome do arquivo
        open OUT, '>', $output_file or die "Não foi possível criar o arquivo '$output_file': $!";
        print OUT "$line\n";  # Escreve o cabeçalho da nova sequência
    } else {
        print OUT "$line\n";  # Escreve a sequência
    }
}

close $fh;
close OUT if $output_file;  # Fecha o último arquivo se existir
