#! /usr/bin/perl
####################
# collapse_snpEff.pl - for a snpEff output collapse snps with multiple lines (transcript) into one. 
# By Hua Li
# 5/2011
####################

#call the script with: perl collapseSnpEff.pl snpEffOutput.txt> output.txt

#for each line in the file
		
$id = '';
$TranscriptAll = '';
$EffectAll = '';
$OldAAn2NewAll = '';
$OldCodon2NewAll = '';
$lastID='';


print ("Chrom", "\t", "Position","\t", "Reference", "\t", "Change", "\t", "ChangeType", "\t", 
       "Quality", "\t", "Coverage", "\t", "GeneID", "\t", "GeneName", "\t", "BioType", "\t", 
       "TranscriptID", "\t", "Effect", "\t", "OldAA2NewAA", "\t", "OldCondon2NewCondon", "\n");

#$lastPosition\t$lastReference\t$lastChange\t$lastChangeType\t$lastHomozygous\t$lastQualityCoverage\t$lastGeneID\t$lastGeneName\t$lastBioType\t$TranscriptAll\t$EffectAll\t$OldAAn2NewAll\t$OldCodon2NewAll\n";
while( <>)
 
{
	if($_ =~ /^#/)
	{
		#print "$_";
		$lastID = ''; 
	}
	else
	{
		chomp;

		($Chrom,$Position,$Reference,$Change,$ChangeType,$Homozygous,$Quality, $Coverage,$Warnings,$GeneID,$GeneName,$BioType,$TrancriptID, $ExonID, $ExonRank, $Effect,  $oldAAn2NewAA,$OldCodon2NewCodon,$CodonNum,$CDS_size, $CodonsAround,$AAsAround, $CustomIntervalID) = split("\t",$_);

		#($chrom,$start,$end,$id,$score,$strand,$pos,$count) = split("\t",$_);
		
		$id = join(';', $Chrom,$Position);

		if($id eq $lastID)
		{
			#keep going
			#$last_count = $last_count + $count;
			$TranscriptAll = join(';',$TranscriptAll, $TrancriptID);
			$EffectAll = join(';', $EffectAll, $Effect);
			$OldAAn2NewAll = join(';', $oldAAn2NewAll, $oldAAn2NewAA);
			$OldCodon2NewAll = join(';', $OldCodon2NewAll, $OldCodon2NewCodon);
			
			#print "$TranscriptAll\n";
		}
		else
		{
		        if ($lastID eq '' and $id ne '')
			{
				$TranscriptAll = $TrancriptID;
				$EffectAll = $Effect;
				$OldAAn2NewAll = $oldAAn2NewAA;
				$OldCodon2NewAll = $OldCodon2NewCodon;	
			}
			
			elsif ( $id ne $lastID)
			{
				print "$lastChrom\t$lastPosition\t$lastReference\t$lastChange\t$lastChangeType\t$lastQuality\t$lastCoverage\t$lastGeneID\t$lastGeneName\t$lastBioType\t$TranscriptAll\t$EffectAll\t$OldAAn2NewAll\t$OldCodon2NewAll\n";

                                $lastID = $id;
				$TranscriptAll = $TrancriptID;
				$EffectAll = $Effect;
				$OldAAn2NewAll = $oldAAn2NewAA;
				$OldCodon2NewAll = $OldCodon2NewCodon;	
			}
		}
		$lastID = $id;
		$lastChrom=$Chrom;
		$lastPosition = $Position;
		$lastReference = $Reference;
		$lastChange = $Change;
		$lastChangeType = $ChangeType;
		#$lastHomozygous = $Homozygous;
		$lastQuality = $Quality;
		$lastCoverage=$Coverage;
		$lastGeneID = $GeneID;
		$lastGeneName = $GeneName;
		$lastBioType = $BioType; 
		
	}
}

print "$lastChrom\t$lastPosition\t$lastReference\t$lastChange\t$lastChangeType\t$lastQuality\t$lastCoverage\t$lastGeneID\t$lastGeneName\t$lastBioType\t$TranscriptAll\t$EffectAll\t$OldAAn2NewAll\t$OldCodon2NewAll\n";


