#!/usr/local/bin/perl
my($stem, $commentchar) = @ARGV;
open STEP, "$stem.step" or die "opening $stem.step";
open MK, ">$stem.mk" or die "opening $stem.mk";
open SCRIPT, ">$stem" or die "opening $stem";
while (<STEP>)
{ if (/$commentchar\s*requires:\s*(\S.*)$/i)
  { print MK "$stem.out : $1\n"; }
  elsif (/$commentchar\s*produces:\s*(\S.*)$/i)
  { print MK "$1 : $stem.out ;\n"; }
  print SCRIPT;
}
