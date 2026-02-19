#!/usr/bin/perl
use strict;
use warnings;

my $input = shift or die "Uso: perl remover_8A.pl <arquivo_multifasta>\n";

open my $IN, '<', $input or die "Não foi possível abrir $input: $!";
open my $OUT_KEEP, '>', 'sem_8A.fasta'    or die "Erro ao criar sem_8A.fasta: $!";
open my $OUT_REMOVED, '>', 'com_8A.fasta' or die "Erro ao criar com_8A.fasta: $!";

# Lê registros FASTA completos, sem perder novoslines
local $/ = "\n>";  

while (my $entry = <$IN>) {
    chomp $entry;
    $entry =~ s/^>//;    # remove '>' se estiver no início

    # Separa header e sequência, preservando quebras de linha originais
    my ($header, @rest) = split /\n/, $entry;
    my $sequence = join("\n", @rest);

    # Reconstrói o registro FASTA exatamente
    my $record = ">$header\n$sequence\n";

    # Filtra pelo padrão 8.A no header
    if ($header =~ /\b8\.A\b/) {
        print $OUT_REMOVED $record;
    } else {
        print $OUT_KEEP    $record;
    }
}

close $IN;
close $OUT_KEEP;
close $OUT_REMOVED;
