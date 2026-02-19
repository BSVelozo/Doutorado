#!/usr/bin/perl
use strict;
use warnings;

# Verifica se os argumentos foram passados corretamente
if (@ARGV != 3) {
    die "Uso: $0 <arquivo_identificadores.3line> <arquivo_multifasta.fasta> <arquivo_resultado.fasta>\n";
}

# Pega os arquivos passados como argumentos
my ($identifiers_file, $fasta_file, $output_file) = @ARGV;

# Função para extrair identificadores e sufixos do arquivo .3line
sub extract_identifiers {
    my ($file) = @_;
    my %identifiers;
    
    open my $fh, '<', $file or die "Não foi possível abrir o arquivo '$file': $!";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^(>\S+)\s+(SP\+TM|TM)$/) {
            my ($id, $suffix) = ($1, $2);
            $identifiers{$id} = $suffix;
        }
    }
    close $fh;
    return %identifiers;
}

# Função para processar o arquivo multifasta e gerar o arquivo de saída
sub process_fasta_file {
    my ($fasta_file, $identifiers_ref, $output_file) = @_;
    my %identifiers = %{$identifiers_ref};
    
    open my $fasta_fh, '<', $fasta_file or die "Não foi possível abrir o arquivo '$fasta_file': $!";
    open my $output_fh, '>', $output_file or die "Não foi possível criar o arquivo '$output_file': $!";
    
    my $write_flag = 0;
    while (my $line = <$fasta_fh>) {
        chomp $line;
        if ($line =~ /^(>\S+)/) {
            my $identifier = $1;
            if (exists $identifiers{$identifier}) {
                print $output_fh "$identifier $identifiers{$identifier}\n";
                $write_flag = 1;
            } else {
                $write_flag = 0;
            }
        } elsif ($write_flag) {
            print $output_fh "$line\n";
        }
    }
    
    close $fasta_fh;
    close $output_fh;
}

# Extrair os identificadores e sufixos
my %identifiers = extract_identifiers($identifiers_file);

# Processar o arquivo multifasta e escrever os resultados
process_fasta_file($fasta_file, \%identifiers, $output_file);

print "Processamento concluído. Resultados salvos em '$output_file'.\n";
       
