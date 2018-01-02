#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Seq;
use File::Copy;

sub help{
	print  "Busqueda de  similitud en bases de datos HMM para identificar familias de proteinas o dominios.
	Opciones:
	-i	Archivo de contigs o proteinas. Default: contigs.fasta
	-d	Base de datos. Default: /data/pfamtigrfam/pfamatigrfam.hmm
	-e	Umbral de e-value permitido. Default: 10
	-h 	Imprime este mensaje de ayuda\n\n"
}

my %options;
getopt('iden', \%options);

#Parametros configurables:	
my %params = (
	i => 'contigs.fasta',				#Sequencias de entrada por defecto
	d => '/data/pfamtigrfam/pfamatigrfam.hmm', 	#Base de datos de Pfam-A +TIGRfam
	e => '10',												#Umbral de e-value
);

#Asigno los parametros que hayan sido definidos externamente
foreach my $option (keys %options) {
	if ($option eq 'h') {&help; exit 1;}
	$params{$option} = $options{$option} if exists $params{$option};
}

my $log = $params{"i"} . ".log";
my $paso_3 =  '';
open(LOG, ">>$log");



if ($params{"i"} =~ /\.faa$/){
	print LOG "Comenzando Pipeline desde paso 3 (usando set de proteinas preexistente).\n";
	print "Iniciando el pipeline desde el paso 3  (usando set de proteínas preexistente).\n";
	$paso_3 = "True";
}
else{
	my $warning = "";
	if (`grep \"Paso 1\" $log -c` == 0) {
		print "Advertencia: el Paso 1 no parece haber sido ejecutado.\n";
		$warning = "True";
	}
	if (`grep \"Paso 2\" $log -c` == 0) {
		print "Advertencia: el Paso 2 no parece haber sido ejecutado.\n";
		$warning = "True";
	}
	#if (`grep \"Paso 3\" $log -c` >= 0) {
	#	print "Advertencia: ya se ejecuto el Paso 3.\n";
	#}

	if ($warning eq "True") {
		print "Ejecutar de todos modos (s/n)?";
		my $seguir = <STDIN>;
		chomp($seguir);
		if ($seguir eq 's'){
			print "";
		}
		else {
			exit 0;
		}
	}
}


#Chequeo la existencia del archivo fasta
my $ls = `ls '$params{"i"}'  2>&1`;
chomp($ls);
unless ($ls eq $params{"i"}){
	die "Error: No se encuentra el archivo $params{\"i\"} \n";
}

print "Busqueda de familias y dominios con HMMER\n";
if ($paso_3 eq "True") {

	my @mock_contig = split("\\.", $params{"i"});
	my $mock_seq = Bio::Seq->new( -seq => 'NNNNNNNNNNNNNNNNNNNNNN', -id  => $mock_contig[0]);
	mkdir $mock_seq->id;
	mkdir "./tmpbia";
	copy($params{"i"},   $mock_seq->id . "/orfs_" . $params{"i"});

	$params{"i"} = $mock_contig[0] .".fna";
	my $mock_genome =  Bio::SeqIO->new( -file => ">$params{\"i\"}", -format => 'fasta');

	$mock_genome->write_seq($mock_seq);
}


my $genomeparser = Bio::SeqIO->new( -file   => $params{"i"}, -format => "fasta");

while (my $seq_record = $genomeparser->next_seq){
	my @hits;
	print "	Procesando " . $seq_record->id . "\n";
	my $hmmerout =  "tmpbia/orfs_" . $seq_record->id . ".hmmer";
	eval{	my $orffile = $seq_record->id . "/orfs_" . $seq_record->id . ".faa";
		`hmmscan -E $params{"e"} --notextw   -o $hmmerout '$params{"d"}' '$orffile' `;

		#Parseo del archivo fasta original
		my $queryparser = Bio::SeqIO->new( -file   => $seq_record->id . "/orfs_" . $seq_record->id . ".faa", -format => "fasta");

		#Parseo del archivo hmmer generado
		my $hmmerparser = Bio::SearchIO->new(-file   => $hmmerout, -format => 'hmmer3');

		while (my $result = $hmmerparser->next_result){				#Obtengo el siguiente resultado de hmmer
			my $queryrec = $queryparser->next_seq;					#Obtengo la siguiente secuencia
			while(my $model = $result->next_model){					#Para cada modelo (hit)
				my $sbjctID  = $model->name;							#Obtengo ID de cada modelo subject
				my $sbjttProd = $model->description;					#Obtengo el producto
				#print $model->significance, "\n";						#E-value global de cada familia								
				while (my $hsp = $model->next_hsp){					#Para cada dominio (hsp) de cada familia
					my $cov = int(100*($hsp->end('query') - $hsp->start('query'))/length($queryrec->seq));	#Calculo la cobertura
					my $ident = ($hsp->frac_identical)*100;			#Porcentaje de identidad
					my $score = $hsp->score;							#Score
					my $evalue = $hsp->evalue;							#E-value								
					my $dfrom =  $hsp->start('query');					#Posicion de inicio y fin del dominio en el query
					my $dto =  $hsp->end('query');
					my $tfrom =  $hsp->start('hit');						#Posicion de inicio y fin del dominio en el modelo
					my $tto =  $hsp->end('hit');
					#Genero un array con los datos del parseo para cada alineamiento
					push (@hits, [$queryrec->id,$queryrec->desc,$sbjctID,$sbjttProd,$dfrom,$dto,$tfrom,$tto,$ident,$cov,$evalue,$score]);
				}	
			}
		}
	};

	my $input_name = "orfs_" . $seq_record->id  . ".hmmer";
	my $hmmeroutfinal = $seq_record->id . "/" . $input_name;

	open(SALIDA_HMMER,">$hmmeroutfinal");

	#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
	print SALIDA_HMMER "ID del query        Descripcion        ID del hit        Producto        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value        score\n";

	my $temp;
	foreach (@hits){
		$temp =  join("\t",@$_[0..11]) . "\n";
		print SALIDA_HMMER  $temp;
	}	


	close(SALIDA_HMMER);

};

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

$year += 1900;
$mon++;

my $hora =  "$mday/$mon/$year $hour:$min:$sec";

print LOG "[ $hora ] Paso 3 HmmerHits corrido sobre archivo $params{\"i\"} usando $params{\"d\"} como base de datos de Modelos de Markov.\n";

close(LOG);
