# change line continuation. If a line has '.' or '>' on the
# fifth position, change into space and at '&' to end of previous line

# quit unless we have an input filename
$num_args = $#ARGV + 1;

if ($num_args != 1) {
    print "\nUsage: comment77to90 <inputfile.f>\n";
    exit;
}

$filename=$ARGV[0];

# ****** count file lines **********

open(filepointer, $filename) or die ("input file not found");

$nlines = 0;
foreach $line (<filepointer>)
{
 $nlines++;
}

#print "FP contains $nlines lines\n";

close(filepointer);

#**** read file and store *********


@text = '';


open(filepointer, $filename) or die ("input file not found");

foreach $line (<filepointer>)
{
    $line =~ s/\R//g;
    @chars = split('',$line);
    push (@text, [@chars])    ;
}


#printf "@text\n";

close(filepointer);


#********* delete old file ********

unlink $filename;

#make a new file with the same name:

unless(open NEWFILE, '>'.$filename) {
    die "nUnable to create $filen";
}


#******** execute corrections ****

@newtext = '';

@previousline = '';
@thisline = '';



    for $i ( 1 .. $#text ) {

        @previousline = @thisline;
        @thisline = '';

for $j ( 0 .. $#{$text[$i]} ) {
            push (@thisline, $text[$i][$j]);
        }    

#        print "previous: @previousline\n";
#        print "this: @thisline\n";

        # make string to write to file
        $string = '';

$nchars = scalar @previousline;
for ($j=1; $j <= $nchars; $j++)
{
    $string = $string.@previousline[$j];
}

#        print("writing: $string");


        if (@thisline[6] eq '>' or @thisline[6] eq '.' or @thisline[6] eq '<' or @thisline[6] eq ',')
        {
            $string = $string.' &';
            @thisline[6] = ' ';

#           print "new previous: $string\n";
#           print "new this: @thisline\n";


        }

        if ($i > 1)   # to prevent first empty line, because no previous line is known, yet
        {
            print NEWFILE "$string\n";
        }

 }
    
    # and finished the last thisline
    $nchars = scalar @thisline;
    $string = '';
 for ($j=1; $j <= $nchars; $j++)
 {
$string = $string.@thisline[$j];
 }

    print NEWFILE "$string\n";
