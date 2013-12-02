#!/usr/bin/perl -w

use DBI;

$SIG{__WARN__} = sub { die @_ };

my $mysql = DBI->connect("DBI:mysql:database=seq_demo;user=seq_demo;password=demo_pass");

$mysql->do(q{LOCK TABLES prot WRITE,
	     annot WRITE,
	     sp WRITE });

my $EL = 125;
my $NA = 123;

my @aatrans = ($EL,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$EL,$NA,$NA,$EL,$NA,$NA,
	       $NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,
	       $NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA, 24,$NA,$NA,$NA,$NA,$NA,
	       $NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,$NA,
	       $NA,  1, 21,  5,  4,  7, 14,  8,  9, 10,$NA, 12, 11, 13,  3,$NA,
	        15,  6,  2, 16, 17,$NA, 20, 18, 23, 19, 22,$NA,$NA,$NA,$NA,$NA,
	       $NA,  1, 21,  5,  4,  7, 14,  8,  9, 10,$NA, 12, 11, 13,  3,$NA,
	        15,  6,  2, 16, 17,$NA, 20, 18, 23, 19, 22,$NA,$NA,$NA,$NA,$NA
	      );

my $ins_prot = $mysql->prepare(q{
    INSERT INTO prot (seq,bin,len) VALUES (?, ?, ?)
    });

my $ins_annot = $mysql->prepare(q{
    INSERT INTO annot (gi, prot_id, db, descr) VALUES (?, ?, ?, ?)
    });

my $ins_sp = $mysql->prepare(q{
    INSERT INTO sp (gi, acc, name) VALUES (?, ?, ?)
    });

use vars qw( $seq $bin $tot_seq $tot_annot $tot_sp );
use vars qw( $gi $prot_id $db $desc $sp_acc $sp_name );
use vars qw( $header $seq @entries );
use vars qw( $gi $db $db_acc $db_name $desc);

$tot_seq = $tot_annot = $tot_sp = 0;

for my $db_file ( @ARGV ) {
    open(DATA, "<$db_file") or die $!;
    local $/ = "\n>";
    while (<DATA>) {
	chomp; # remove trailing "\n>" record header
	($header, $seq) = $_ =~ m/^>?  # record separator (first entry)
	    ( [^\n]* ) \n  # header line
		(     .* )     # the sequence
		    /osx; # optimize, multiline, commented
	
	$seq =~ s/\W|\d//sg;
	$bin = pack('C*', map { $aatrans[unpack('C', $_)] } split(//, $seq));
	$ins_prot->execute($seq,$bin,length($seq));
	$prot_id = $ins_prot->{mysql_insertid};

	$tot_seq++;

#	print STDERR "Inserted $prot_id: ". length($seq)."\n";

	@entries = split(/\001/, $header);

	for ( @entries ) {
	    ($gi,$db,$db_acc,$db_name,$desc)= 
		$_ =~ /^gi\|(\d+)\|([a-z]+)\|(\S*)\|(\S*) (.*)$/o;
#	    print "$prot_id: $gi\t$db\t$db_acc\t$desc\n";
	    $ins_annot->execute($gi,$prot_id,$db,$desc);

	    $tot_annot++;

	    if ($db eq "sp") {
		$ins_sp->execute($gi,$db_acc,$db_name);
		$tot_sp++;
	    }
	}
    }
    close(DATA);
}

print "Inserted $tot_seq sequences; $tot_annot annotations; $tot_sp swissprot\n";



