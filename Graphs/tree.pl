#!usr/local/bin/perl

open(DOT,">graph.dot");
open(LU,"+>>lookup.temp");


$count =0;
$FoundFile = 0;
$ext="f90";
my @files0 = glob("*.$ext");
my @files1 = glob("*/*.$ext");

@files = (@files0,@files1);

print DOT "digraph G { \n
	nodesep=.05; \n
	rankdir=LR; \n
	node [shape=record,width=.1,heigth=.1];\n";

#print FH "nodeL [ label = \"{A-M}\"];\n
#	nodeR [ label = \"{N-Z}\"];\n";
$filesparsed = 0;
foreach $file (@files) 
{
	open(Search,"<$file");
	while(<Search>)
		{
		if($_ =~ "call")
			{
			$StoreName = $_;
			($Store) = ($StoreName =~ /^\s*call ([\d\w]+)/);
			if(length($1)<2){next;}
			#$Store = $Store."$ext";
			#($Store) = ($Store =~ /\/(.*)$ext/);
			#$Store = $Store."$ext";
			#substr($Store,0,23)="";
			$filename = $file;
			$filename=~s/(?:.*\/)?([^\/]+)\.$ext$/$1/;
			#substr($filename,0,4)="";
			print "$filename,$Store\n";
			print LU "$filename,$Store\n";
			open(UsedRead, "<used.temp");
			$FoundFile = 0;
			$FoundStore =0;
			while(<UsedRead>)
				{
				chomp;
				($UsedFlag, $NodeStuff) = split(/,/);
				if($UsedFlag eq $filename) {$FoundFile=1;}
				if($UsedFlag eq $Store) {$FoundStore=1;}
				}
			close(UsedRead);
			if($FoundFile == 0)
				{
		                print DOT "node$count [label = \"$filename\"];\n";
				open(UsedWrite, ">>used.temp");
				print UsedWrite "$filename,node$count\n";
				
				$count++;
				close(UsedWrite)
				}
			if($FoundStore == 0)
				{
				print DOT "node$count [label = \"$Store\"];\n";
				open(UsedWrite, ">>used.temp");
				print UsedWrite "$Store,node$count\n";
				$count++;
				close(UsedWrite)
				}
			

			}
		}

		$filesparsed++;
	close(Search);
}

#Now make dot graph
#
close(LU);
open(MakeDot,"<lookup.temp");
while(<MakeDot>)
	{
	chomp;
	($first, $second) = split(/,/);
	open(GetNodeNums,"<used.temp");
	while($UsedGet = <GetNodeNums>)
		{
		chomp $UsedGet;
		($NodeName, $NodeNum) = split(/,/,$UsedGet);
		if($NodeName eq $first)
			{
			$DotLeft = $NodeNum;
			}
		if($NodeName eq $second)
			{
			$DotRight = $NodeNum;
			}
		}
	close(GetNodeNums);
	print DOT "$DotLeft -> $DotRight;\n";

	}
		

#print DOT "nodetotal [label = \"$count calls; $filesparsed files\",color=red];\n";
print DOT "}\n";
close(DOT);
close(MakeDot);
close(Used);
unlink("lookup.temp");
unlink("used.temp");
