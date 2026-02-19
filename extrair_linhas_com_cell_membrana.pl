#!/usr/bin/perl
use strict;
use warnings;

# Verifica se o arquivo de entrada foi passado como argumento
if (@ARGV != 1) {
    die "Uso: $0 arquivo.csv\n";
}

my $input_file = $ARGV[0];
my $output_file = "linhas_membrana.csv";

# Abre o arquivo de entrada para leitura
open(my $in_fh, '<', $input_file) or die "Não foi possível abrir o arquivo '$input_file': $!";

# Abre o arquivo de saída para escrita
open(my $out_fh, '>', $output_file) or die "Não foi possível criar o arquivo '$output_file': $!";

# Lê a primeira linha (cabeçalho)
my $header = <$in_fh>;
chomp $header;

# Divide o cabeçalho pelos separadores de vírgula
my @headers = split(",", $header);

# Verifica se a segunda coluna é "Localizations"
if ($headers[1] ne "Localizations") {
    die "A segunda coluna não é 'Localizations'. Verifique o formato do arquivo.\n";
}

# Escreve o cabeçalho no arquivo de saída
print $out_fh "$header\n";

# Percorre cada linha do arquivo de entrada
while (my $line = <$in_fh>) {
    chomp $line;

    # Divide a linha pelos separadores de vírgula
    my @columns = split(",", $line);

    # Verifica se a célula da coluna "Localizations" contém "Cell membrane"
    if (defined $columns[1] && $columns[1] =~ /Cell membrane/i) {
        print $out_fh "$line\n";  # Escreve a linha no arquivo de saída
    }
}

# Fecha os arquivos
close($in_fh);
close($out_fh);

print "Linhas com 'Cell membrane' na coluna 'Localizations' foram salvas no arquivo '$output_file'.\n";
