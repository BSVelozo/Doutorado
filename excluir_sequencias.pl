#!/usr/bin/perl
use strict;
use warnings;

#Esse script compara uma lista de ids com um arquivo multifasta e exlui as sequencias que possuem os ids daquela lista

# Função para carregar identificadores da lista e adicionar '>'
sub carregar_identificadores {
    my ($arquivo_lista) = @_;
    my %identificadores;

    open(my $fh, '<', $arquivo_lista) or die "Não foi possível abrir o arquivo '$arquivo_lista' $!";
    while (my $linha = <$fh>) {
        chomp($linha);
        $identificadores{">" . $linha} = 1;  # Adiciona '>' ao início de cada identificador
    }
    close($fh);

    return %identificadores;
}

# Função para excluir sequências cujos cabeçalhos estão na lista de identificadores
sub excluir_sequencias {
    my ($arquivo_fasta, $identificadores_ref, $arquivo_saida) = @_;
    my %identificadores = %{$identificadores_ref};

    open(my $fasta_fh, '<', $arquivo_fasta) or die "Não foi possível abrir o arquivo '$arquivo_fasta' $!";
    open(my $saida_fh, '>', $arquivo_saida) or die "Não foi possível criar o arquivo '$arquivo_saida' $!";

    my $escrever = 1;
    while (my $linha = <$fasta_fh>) {
        if ($linha =~ /^>/) {
            # Verifica se o cabeçalho está na lista de identificadores
            my $cabecalho = $linha;
            chomp($cabecalho);
            $escrever = !exists $identificadores{$cabecalho};
        }
        print $saida_fh $linha if $escrever;
    }

    close($fasta_fh);
    close($saida_fh);
}

# Arquivo com a lista de identificadores
my $lista_identificadores_file = 'Id_warning.txt';
# Arquivo FASTA de entrada
my $fasta_file = 'gencode.fasta';
# Arquivo FASTA de saída com as sequências filtradas
my $output_file = 'sequencias_filtradas.fasta';

# Carregar identificadores da lista
my %identificadores = carregar_identificadores($lista_identificadores_file);

# Excluir sequências cujos cabeçalhos estão na lista
excluir_sequencias($fasta_file, \%identificadores, $output_file);

print "Processo concluído. Sequências filtradas estão no arquivo: $output_file\n";
