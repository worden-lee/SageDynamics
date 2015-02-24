#!/usr/local/bin/perl
my($stem, $commentmark) = @ARGV;
open STEP, "$stem.step" or die "opening $stem.step";
open MK, ">$stem.mk" or die "opening $stem.mk";
open SCRIPT, ">$stem" or die "opening $stem";
print MK "# helper makefile automatically generated from $stem.step\n";
my $products = '';
while (<STEP>)
{ if (/$commentmark\s*requires:\s*(\S.*)$/i) {
    print MK "$stem.out $stem.tried : $1\n";
  }
  elsif (/$commentmark\s*produces:\s*(\S.*)$/i) {
    #print MK "$1 : $stem.tried ;\n";
    print MK "$1 : $stem.out ;\n";
    $products .= " $1";
  }
  print SCRIPT;
}
if ($products) {
  print MK "$stem.out $stem.tried : STEP_PRODUCTS=$products\n";
}
