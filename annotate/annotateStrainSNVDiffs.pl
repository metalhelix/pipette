#! /usr/bin/env perl

# by Aaron Noll originally
# edited by Madelaine Gogol 3/2011

use strict;
if (@ARGV < 2)
{
	die "usage: $0 vcfFile sliceAdaptorName dbName port host";
}

my ($vcfFile,$sliceAdaptorName,$dbName,$port,$host) = @ARGV;

use lib "/home/solexa/ensembl/src/ensembl_57/modules";
use Bio::EnsEMBL::Registry;
use DBI;
use Bio::PrimarySeq;
use Data::Dumper;
use POSIX;

my $impactReportHeader =  "Gene,Transcript(pos;strand),Impact,GO,InterPro\n";

open(vcfFile,$vcfFile);
my %nucInfo = ("A","","C","","G","","T","");
sub getMutPos($$);
my $goDB = "go_201001";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					-host => $host,
					-DBNAME => $dbName,		
					-user => 'anonymous',
					-port => $port 
) or die "Could not connect to database";

my $goh = DBI->connect('DBI:mysql:go_201001:mysql-dev', 'anonymous') or die "Could not connect to GO database";;

print STDERR "$sliceAdaptorName\n";

my $sa = $db->get_SliceAdaptor($sliceAdaptorName,'Core','Slice');
my $slice;

while(<vcfFile>)
{
	my @goData = ();
	my @interProData = ();
	my $sep = "; ";	
	chomp();

	if($_ =~ /^##/)
	{
	}
	elsif($_ =~ /^#/)
	{
		my $header = $_;	
		$header =~ s/\t/,/g;
		print "chr,start,end,ref,alt,qual," . $impactReportHeader; 
	}
	else
	{
		#adding a tab on end of line for later
		#my $line = $_ . "\t";
		$_ =~ s/"//g;

		my @line = split("\t",$_);
		my $chr = $line[0];
		my $chrPos = $line[1];
		my $ref = $line[3];
		my @mutAlleles = split(",",$line[4]);
		my $qual = $line[5];
		my $reflen = length($ref);
		my $end = $chrPos + $reflen - 1;

		print STDERR "chr:$chr chrPos:$chrPos end:$end\n";

		$slice = $sa->fetch_by_region('chromosome', $chr, $chrPos, $end);

		#get the sequence from the slice
		my $seq = $slice->seq();
		print STDERR "seq:$seq\n";

		#get some features from the slice
		my %nucs = ();
		foreach my $mutAllele(@mutAlleles)
		{
			if ($mutAllele eq $ref)
			{
				next;
			}
			$nucs{$mutAllele} = "";
		}
		if (scalar(keys %nucs) == 0)
		{
			next;
		}
		print STDERR "going to print Dumper for nucs (mutAlleles)";
		print STDERR Dumper(%nucs);
		my @transcriptColumn = ();
		my @geneColumn = ();
		my @transcripts = @{$slice->get_all_Transcripts()};
		my @exons = @{$slice->get_all_Exons()};
		my $geneAdaptor; 
		my $gene;
		if ((scalar @transcripts) > 1)
		{
			$sep = '; ';
		}
		foreach my $transcript (@transcripts)  #this runs for each transcript.
		{
			my $geneAdaptor = $db->get_GeneAdaptor;
			my $gene = $geneAdaptor->fetch_by_transcript_stable_id($transcript->stable_id);	
			push (@geneColumn,$gene->external_name);
			#push (@geneColumn,$gene->display_id); #this is FBID
			my $strand = $transcript->strand();
			print STDERR "strand:$strand\n";
			#only spliced exons - nothing more
			my $tseq = $transcript->translateable_seq();	
			#protein before snp
			my $eID = $transcript->external_name;	
			print STDERR "transcript:$eID\n";	
			print STDERR "Transcript is on slice: ", $transcript->slice->name, "\n";
			print STDERR "Transcript cdsCoords: ", $transcript->start, ' to ', $transcript->end , "\n";
			print STDERR "Transcript seqCoords: ", $transcript->seq_region_start , ' to ', $transcript->seq_region_end, "\n";
			my $len = $transcript->length;
			my $trmaper = $transcript->get_TranscriptMapper();

			#start at genomic coordinate of 1, get only 1, and base
			#has 5' UTR and 3' UTR
			my @cdnaCoords = $trmaper->genomic2cdna(1,1,$strand);
			#only exons
			my @cdsCoords = $trmaper->genomic2cds(1,1,$strand);
			print STDERR "transcript cdsCoord " . Dumper(@cdsCoords);
			print STDERR "transcript cdnaCoord " . Dumper(@cdnaCoords);
			my $snpPos = 0;
			if ($transcript->biotype eq 'protein_coding')
			{
				print STDERR "protein coding \n";
				if ($cdsCoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'))
				{
					my $snpPos = $cdsCoords[0]->start;
					printf STDERR "%s %s -> %s non-gap\n", $transcript->stable_id, $chrPos, $cdsCoords[0]->start, $cdsCoords[0]->strand;
					print STDERR "Transcript length: ", $len,"\n";
					print STDERR "transcript sequence before SNP is ",$tseq . "\n";	
					print STDERR "transcript display ID " . $transcript->translation()->display_id() . "\n";	
					my $oprot = $transcript->translation()->seq();	
					#VALUE COMES OUT TO BE = ($transcript->length - $transcript->end) 
					#NB subtract 1 later to make 0 based when making substitution function
					print STDERR "snpPos is $snpPos \n";
					push (@transcriptColumn,"$eID($snpPos;$strand)");
					#look at each strain's SNPs
					foreach my $nuc (sort keys %nucs) 
					{
						#position for SNP zero indiced for substr()... build new seq with new nucleotide
						if ($strand == 1)
						{
							print STDERR "strand: $strand... if given ref then trancript should change at $snpPos from $ref to $nuc\n";
						}
						else
						{
							my $tmpRef = $ref;
							my $tmpNuc = $nuc;
							$tmpNuc =~ tr/ATGCatgc/TACGtacg/;
							$tmpRef =~ tr/ATGCatgc/TACGtacg/;
							print STDERR "strand: $strand... if given ref then trancript should change at $snpPos from $tmpRef to $tmpNuc\n";
						}
						my $nseqstr = $tseq;	
						#make substitution in nucleotide sequence
						if ($strand < 0)
						{
							#temporarily change for negative strand but change back later for display purposes
							$nuc =~ tr/ACGT/TGCA/;	
							substr($nseqstr,($snpPos-1),1,$nuc);
							$nuc =~ tr/ACGT/TGCA/;	
						}
						else
						{
							substr($nseqstr,($snpPos-1),1,$nuc);
						}
						my $nseq = Bio::PrimarySeq->new(-seq => $nseqstr);
						print STDERR "strand: $strand... trancript changed from " . substr($tseq,($snpPos-1),1) . " to ". substr($nseqstr,($snpPos-1),1) . "\n";
						print STDERR "strand: $strand... transcript sequence after SNP is ",$nseqstr . "\n";	
						$len = length($oprot);
						print STDERR "Protein length: ", $len,"\n";
						#get new protein sequence 
						my $nprot = $nseq->translate()->seq();
						print STDERR "protein sequence before SNP $oprot\n";	
						print STDERR "protein sequence after SNP $nprot\n";	
						if (substr($nprot,-1) eq '*')
						{	
							#remove stop codon since ensembl doesn't include
							chop($nprot);
						}
						else
						{
							#lost stop codon!!
							$nucs{$nuc} = "LS(*".length($nprot). substr($nprot,-1) .")";	
							next;
						}
						if ($oprot ne $nprot)
						{
							my $stopPos = index($nprot,'*');
							if ($stopPos >0)
							{
								$nucs{$nuc} = "nonsense($stopPos:".length($nprot).")";
							}
							else
							{
								#GO ANNOTATION
								foreach my $goTerm (@{$transcript->get_all_DBLinks('GO')})
								{
									print STDERR "go:". $goTerm->display_id . "\n";
								}
								my @dbEntries = @{$transcript->translation()->get_all_DBEntries()};
								foreach my $dbEntry (@dbEntries)
								{
									my $primary_id = $dbEntry->primary_id;
									if ($dbEntry->isa('Bio::EnsEMBL::GoXref'))
									{
										print STDERR "FOUND GOADAPTOR\n\n";
										print STDERR "primary ID $primary_id\n";
										my $go_query = $goh->prepare("select acc,name from $goDB.term WHERE acc = ?");
										$go_query->execute($primary_id);	
										while (my $result = $go_query->fetchrow_hashref())
										{
											push(@goData,$result->{"acc"} . "[" . $result->{"name"} . "]");
										}
									}
								}	
								#interpro annotation
								my @domainFeatures = @{$transcript->translation()->get_all_DomainFeatures()};
								foreach my $domainFeature (@domainFeatures)
								{
									if ($domainFeature->interpro_ac() ne "")
									{
										push(@interProData,$domainFeature->interpro_ac() . "[".$domainFeature->idesc()."]");
									}
								}
								$len = length($oprot);
								print STDERR "Protein length: ", $len,"\n";
								my $snpProteinLoc = getMutPos($oprot,$nprot);
								my @cdsCoords = $trmaper->genomic2pep(1,$len,$strand);
								print STDERR "SNP protein change is " ,substr($oprot,$snpProteinLoc,ceil($reflen/3)) , " to ", substr($nprot,$snpProteinLoc,ceil($reflen/3)),"\n";	
								$eID = $gene->stable_id;
								print STDERR $eID . "\n";	
								my $sStr = substr($oprot,$snpProteinLoc,1) . ($snpProteinLoc+1) . substr($nprot,$snpProteinLoc,1);
								$nucs{$nuc} = $sStr; 
								print STDERR "SNP protein location is " , ($snpProteinLoc+1) ,"\n";
							}
						}
						else
						{
							$nucs{$nuc} = "synonymous";
						}
					}
					foreach my $nuc (sort keys %nucs)
					{
						$nucInfo{$nuc} = $nucInfo{$nuc}.=$nucs{$nuc}.$sep;
					}
				}
				elsif($cdnaCoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'))
				{
					$snpPos = $cdnaCoords[0]->start;
					push (@transcriptColumn,"$eID: ". "UTR". "($snpPos;$strand)");
					printf STDERR "%s %s -> %s non-gap\n", $transcript->stable_id, $chrPos, $cdnaCoords[0]->start, $cdnaCoords[0]->strand;
				}
				elsif($cdnaCoords[0]->isa('Bio::EnsEMBL::Mapper::Gap'))
				{
					my ($intron,$msg,$exon);	
					my @introns = @{$transcript->get_all_Introns};
					foreach my $tmp (@introns)
					{
						$intron = $tmp;
						print STDERR $intron->seq_region_start . " start \n";
						print STDERR $intron->seq_region_end . " end \n";
						if ($chrPos == $intron->seq_region_start or $end == $intron->seq_region_start)
						{
							$exon = $intron->prev_Exon();
							if ($intron->strand < 1)
							{
								$exon = $intron->next_Exon();
							}
							$msg = "SJ;Exon:" . $exon->display_id;
							last;
						}
						elsif ($chrPos == $intron->seq_region_end or $end == $intron->seq_region_end)
						{
							$exon = $intron->next_Exon();
							if ($intron->strand < 1)
							{
								$exon = $intron->prev_Exon();
							}
							$msg = "SJ;Exon:" . $exon->display_id;
							last;
						}
						elsif ($chrPos > $intron->seq_region_start && $end < $intron->seq_region_end)
						{
							my $nextExonDist = abs($end - $intron->next_Exon()->seq_region_start);
							my $prevExonDist = abs($chrPos - $intron->prev_Exon()->seq_region_end);
							if ($intron->strand < 1)
							{
								$nextExonDist = abs($chrPos - $intron->next_Exon()->seq_region_end);
								$prevExonDist = abs($end - $intron->prev_Exon()->seq_region_start);
							}
							$msg = "E: " . $intron->prev_Exon()->display_id . " $prevExonDist bps ";
							$msg .= "E: " . $intron->next_Exon()->display_id . " $nextExonDist bps";
							if ($nextExonDist == 20)
							{
								#a branch point is 20bps away from the AG side of the intron (aka the end of the intron)
								my $tmpSeq = $intron->seq();
								my $bpToAGSeq = substr($tmpSeq,length($tmpSeq)-20);
								$msg .= " BP: $bpToAGSeq";
							}	
							last;
						}
					}
					$snpPos = $cdnaCoords[0]->start;
					push (@transcriptColumn,"$eID: ". "Intronic". "($strand) $msg");
					printf STDERR "%s %s -> %s gap\n", $transcript->stable_id, $chrPos, $cdnaCoords[0]->start;
					print STDERR "intron sequence:\n" . $intron->seq() . "\n";
				}
				else 
				{
					push (@transcriptColumn,"$eID: ". "Unknown ");
				}
			}
			else
			{
				push (@transcriptColumn,"$eID: ". $transcript->biotype. "($snpPos;$strand)");
			}
			if (@goData >0)	
			{
				push (@goData,$sep);
			}
			if (@interProData>0)	
			{
				push (@interProData,$sep);
			}
		}
		#NO MORE TRANSCRIPTS

		my @temp = unique(@geneColumn);
		print $chr,",",$chrPos,",",$end,",",$ref,",",join(";",@mutAlleles),",",$qual,",". join($sep,@temp);
		print ",". join($sep,@transcriptColumn);
		print ",";
		foreach my $nuc (sort keys %nucInfo)
		{
			print STDERR "my key is $nuc\n value is ", $nucInfo{$nuc} . "\n";
			print $nucInfo{$nuc};
		}
		print ",\"" . join(",",@goData) ."\"";
		print ",\"" . join(",",@interProData) ."\"";
		print "\n";
		#CLEAR OUT HASH FOR NEXT LINE IN INPUT FILE
		%nucInfo = ("A","","C","","G","","T","");
	}
}

sub getMutPos($$)
{
	my @op = split("",shift(@_));
	my @np = split("",shift(@_));
	#doing this brute force way for now need to find more elegant way
	if ($#op == $#np)
	{
		for(my $i =0;$i<@op;$i++)
		{
			if ($op[$i] ne $np[$i]){return $i;}
		}
	}
	else
	{
		print STDERR "ERROR old and new proteins are different lengths!!";
	}
}

sub unique() #return unique elements of an array
{
	my @input = @_;
	my %seen = ();
	my @uniq = ();
	foreach my $item (@input)
	{
		unless ($seen{$item}) 
		{
			# if we get here, we have not seen it before
			$seen{$item} = 1;         
			push(@uniq, $item);     
		} 
	}
	return(@uniq);
}
