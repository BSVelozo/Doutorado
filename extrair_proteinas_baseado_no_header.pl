#!/usr/bin/perl
use strict;
use warnings;

# Verifica se os arquivos foram passados como argumentos
if (@ARGV != 2) {
    die "Uso: $0 arquivo.csv arquivo.fasta\n";
}

my ($csv_file, $fasta_file) = @ARGV;
my $output_file = "proteinas_extraidas.fasta";

# Abre o arquivo CSV para leitura
open(my $csv_fh, '<', $csv_file) or die "Não foi possível abrir o arquivo CSV '$csv_file': $!";

# Ignora a primeira linha (título da coluna)
my $header_line = <$csv_fh>;

# Lê todos os identificadores de proteínas do CSV e armazena em um hash
my %protein_ids;
while (my $line = <$csv_fh>) {
    chomp $line;
    # Supondo que cada linha tem um identificador na primeira coluna, sem delimitadores específicos
    my @columns = split(/\s+/, $line);  # Divide por espaços ou tabulações
    $protein_ids{$columns[0]} = 1;      # Armazena o identificador da primeira coluna
}
close($csv_fh);

# Abre o arquivo FASTA para leitura
open(my $fasta_fh, '<', $fasta_file) or die "Não foi possível abrir o arquivo FASTA '$fasta_file': $!";

# Abre o arquivo de saída para escrita
open(my $out_fh, '>', $output_file) or die "Não foi possível criar o arquivo de saída '$output_file': $!";

# Variáveis para armazenar temporariamente o header e a sequência da proteína
my $header = '';
my $sequence = '';

# Percorre o arquivo FASTA
while (my $line = <$fasta_fh>) {
    chomp $line;
    
    if ($line =~ /^>(\S+)/) {  # Se a linha começa com '>', é um header
        # Se já temos um header e uma sequência, verificamos o match
        my $identifier = $1;  # Extrai o identificador após o '>'

        if ($header && exists $protein_ids{$header}) {
            # Escreve o header e a sequência no arquivo de saída
            print $out_fh ">$header\n$sequence\n";
        }
        
        # Atualiza o header com o novo identificador
        $header = $identifier;
        $sequence = '';  # Reseta a sequência para a nova entrada
    } else {
        # Concatena a sequência (pode ser dividida em várias linhas no arquivo FASTA)
        $sequence .= $line;
    }
}

# Verifica a última sequência após o loop
if ($header && exists $protein_ids{$header}) {
    print $out_fh ">$header\n$sequence\n";
}

# Fecha os arquivos
close($fasta_fh);
close($out_fh);

print "Proteínas extraídas foram salvas no arquivo '$output_file'.\n";
