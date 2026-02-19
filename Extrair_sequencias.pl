#!/usr/bin/perl
use strict;
use warnings;

# Verificar se os arquivos foram passados como argumento
if (@ARGV != 3) {
    die "Uso: $0 <headers.txt> <multifasta.fasta> <saida.fasta>\n";
}

my ($headers_file, $multifasta_file, $output_file) = @ARGV;

# Ler os headers e armazená-los em um hash (removendo o símbolo '>')
my %headers;
open(my $fh_headers, '<', $headers_file) or die "Não foi possível abrir o arquivo '$headers_file': $!";
while (my $line = <$fh_headers>) {
    chomp $line;
    $line =~ s/^>//;  # Remove o '>' do início da linha
    $headers{$line} = 1 if $line ne '';  # Adicionar somente headers não vazios
}
close($fh_headers);

# Abrir o arquivo multifasta e o arquivo de saída
open(my $fh_multifasta, '<', $multifasta_file) or die "Não foi possível abrir o arquivo '$multifasta_file': $!";
open(my $fh_out, '>', $output_file) or die "Não foi possível criar o arquivo de saída '$output_file': $!";

my $current_header = '';
my $current_sequence = '';

while (my $line = <$fh_multifasta>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
        # Se já temos um header e uma sequência, escrevemos no arquivo de saída
        if ($current_header && exists $headers{$current_header}) {
            print $fh_out ">$current_header\n$current_sequence\n";
        }
        # Atualizar o header atual
        $current_header = $1;
        $current_sequence = '';  # Reiniciar a sequência
    } else {
        # Adicionar a linha à sequência atual
        $current_sequence .= $line;
    }
}

# Verificar se a última sequência deve ser escrita
if ($current_header && exists $headers{$current_header}) {
    print $fh_out ">$current_header\n$current_sequence\n";
}

# Fechar os arquivos
close($fh_multifasta);
close($fh_out);

print "Extração completa! As sequências foram salvas em '$output_file'.\n";
