#!/usr/bin/perl
#parse GE pulse sequence time history file (.tgt.his) and create a matlab script for display

use Data::Dumper qw(Dumper);
use File::Basename;

$fname = basename(shift,'.tgt.his');
$outfile = shift || 'plotpulse.m';   

$fn = $fname.".tgt.his";
open(DATA, "<".$fn);



while(<DATA>) {
     print $_;
     if (/^\s*$/) {
     }
  elsif (/^\*+\s*$/) {
     print "Probably end of Board\n";
     }
  elsif (/^\*Board:\s*(\w+)\s*$/) {
     $board = $1;
     print "Board $board\n";
#     last if $board eq "SSP";
     }
  elsif (/^\s*\*\* Sequence \((\d+)\) \*\*\s*$/) {
     $seq = $1;
     $mc{$board}{$seq} = [];
     print "Sequence $1\n";
     }
  elsif (/^\s*Pulse:\s+(\w+)\s+Type:(\w+)\s+(.*)$/) {
     $pulse = $1;
     $type  = $2;
     if ($board ne "SSP" && $type eq 'Const') {
        if ($3 =~ /^\s*Amp:(-?\d+)\s*$/) {
           $ampS   = $1; 
	   $ampE   = $1;
           }
        else { die "Const line invalid: $_\n"; }
        }
     elsif ($type eq 'Ramp') {
        if ($3 =~ /^\s*Start Amp:(-?\d+)\s+End Amp:(-?\d+)\s*$/) {
	   $ampS   = $1;
	   $ampE   = $2;
           }
        else { die "Ramp line invalid: $_\n"; }
	}
     elsif ($type eq 'Sinc') {
        if ($3 =~ /^\s*Amp:(-?\d+)\s+Nsinc Cycl:(\d+\.\d+)\s+Alpha:(\d+\.\d+)\s*$/) {
           $ampS    = $1;
           $ampE    = $1;
           $ncycles = $2;
           $alpha   = $3;
           }
        else { die "Sinc line invalid: $_\n"; }
        }
     elsif ($type eq 'External') {
        print "External???: $_\n";
        }
     elsif ($type eq 'Reserve') {
	print "Reserve????: $_\n";
	}
     elsif ($type eq 'Bits' || $type eq 'Const') {
        print "Bits???: $_\n";
	chomp;
	my @info = $_;
	$_ = <DATA>; chomp; push(@info, $_);
	if (/^\s+Instr:\s*Start:\s*(\d+)\s+End:\s*(\d+)\s+Dur:\s*(\d+)\s*.*$/) {
	   $start = $1."/1000";
	   $end   = $2."/1000";
	   $dur   = $3."/1000" ;
           }
        while(($_=<DATA>) !~ /^\s*$/) { chomp; push(@info, $_); }
	$y = ($type eq 'Bits') ? 2200 : 2400;
	push(@{$mc{$board}{$seq}}, "s={'".join("'; ...\n'",@info)."'};\n"
		               ."line('XData',[$start,$end],'YData',[".($y-500).",$y]"
			            .",'Color',".($type eq 'Bits' ? "'r'" : "'g'").",'LineWidth',2"
				    .",'Tag','$pulse','ButtonDownFcn',\@(h,e) BDFssp(h,e)"
				    .",'UserData',struct('text',s)"
				    .");");
        next;
        }
     else { die "Unknown type($type): $_\n"; }

     $_ = <DATA>;
     if (/^\s+Instr:\s+Start:(\d+)\s+End:(\d+)\s+Dur:(\d+)\s+Per:\s*(\d+)\s+Amp:\s*(-?\d+)\s*$/) {
	$start = $1."/1000";
	$end   = $2."/1000";
	$dur   = $3."/1000";
	$per   = $4."/1000";
	$amp   = $5;
	}
     else {
	die "Second line of pulse not correct\n";
	}
     $_ = <DATA>;
     if (/^\s+Area:(-?\d+\.\d+)\s+Abs Area:(\d+\.\d+)\s+$/) {
	$area = $1;
	$absarea = $2;
        }
     else {
	die "Third line of pulse not correct\n";
        }
     print "Pulse $pulse $type amp($ampS $ampE $amp) time($start $end $dur $per) \n";
     if ($type eq "Sinc") {
        push(@{$mc{$board}{$seq}}, "line([$start:0.01:$end],$ampS.0*$amp.0*sinc((([$start:0.01:$end] - $start)./($end-$start)-0.5).*4.*$ncycles)"
                                .",'Color','r','LineWidth',2"
                                .",'Tag','$pulse','ButtonDownFcn',\@(h,e) buttondown(h,e)"
                                .",'UserData',struct('Amp',$amp,'AmpS',$ampS,'AmpE',$ampE)"
                                      .");");
        }
     else {
     push(@{$mc{$board}{$seq}}, "line('XData',[$start $end],'YData',[$ampS $ampE].*$amp"
                                .",'Color','r','LineWidth',2"
                                .",'Tag','$pulse','ButtonDownFcn',\@(h,e) buttondown(h,e)"
                                .",'UserData',struct('Amp',$amp,'AmpS',$ampS,'AmpE',$ampE)"
                                .");");
        }
     if ($seq == 1) {
        $xmax = $end if $end > $xmax;
        }
     }
  else {
     print "Invalid line:'$_'\n";
     }
  }

#print Dumper \%mc;

open(my $mf, ">".$outfile) or die "Could not open plotpulse.m\n";
print $mf "function plotpulse()\n";
print $mf "figure;\n";

print $mf "if exist('".$fname.".tgt.pd')\n";
print $mf "    plotPdFile('".$fname.".tgt.pd')\n";
print $mf "    hold on;\n";
print $mf "end\n";

print $mf "subplot(5,1,1);\n";
print $mf "axis([0 $xmax -32768*32768 32768*32768]);\n";
print $mf "line([0 $xmax], [0 0], 'Color', 'b');\n";
print $mf "ylabel('XGRAD');\n";
print $mf join("\n",@{$mc{"XGRAD"}{1}});

print $mf "subplot(5,1,2);\n";
print $mf "axis([0 $xmax -32768*32768 32768*32768]);\n";
print $mf "line([0 $xmax], [0 0], 'Color', 'b');\n";
print $mf "ylabel('YGRAD');\n";
print $mf join("\n",@{$mc{"YGRAD"}{1}});

print $mf "subplot(5,1,3);\n";
print $mf "axis([0 $xmax -32768*32768 32768*32768]);\n";
print $mf "line([0 $xmax], [0 0], 'Color', 'b');\n";
print $mf "ylabel('ZGRAD');\n";
print $mf join("\n",@{$mc{"ZGRAD"}{1}});

print $mf "subplot(5,1,4);\n";
print $mf "axis([0 $xmax -32768*32768 32768*32768]);\n";
print $mf "line([0 $xmax], [0 0], 'Color', 'b');\n";
print $mf "ylabel('RHO1');\n";
print $mf join("\n",@{$mc{"RHO1"}{1}});

print $mf "subplot(5,1,5);\n";
print $mf "axis([0 $xmax 0 2700]);\n";
print $mf "line([0 $xmax], [0 0], 'Color', 'b');\n";
print $mf "ylabel('SSP');\n";
print $mf join("\n",@{$mc{"SSP"}{1}});

print $mf "\n\n";
print $mf "xlh=uicontrol('Style','slider','Min',0,'Max',$xmax,'Value',0,'Units','Normalized','Position',[0.05 0.025 0.2 0.05],'Callback',\@(h,e) xlimCB(h,e));\n";
print $mf "xuh=uicontrol('Style','slider','Min',0,'Max',$xmax,'Value',$xmax,'Units','Normalized','Position',[0.7 0.025 0.2 0.05],'Callback',\@(h,e) xlimCB(h,e));\n";
print $mf "xedt=uicontrol('Style','text','Units','Normalized','Position',[0.35,0.94,0.3,0.04],'String','Pulse info','BackgroundColor',[1,1,1],'FontSize',10);\n";

print $mf "linkaxes(findobj(gcf, 'type', 'axes'),'x');\n";


print $mf "ax = findobj(gcf,'type','axes');\n";
print $mf "for k=1:numel(ax)\n";
print $mf "    addPulseList(ax(k))\n";
print $mf "end\n";

print $mf "function xlimCB(h,e)\nlv=get(xlh,'Value');\nhv=get(xuh,'Value');\nfor sp=1:5\nsubplot(5,1,sp);\nxlim([lv hv]);\nend\nend\n";
print $mf "\n\n";
print $mf "function buttondown(h,e)\n";
print $mf "u=get(h,'UserData');\n";
print $mf "str = sprintf('%s Amp=%d\\n',get(h,'Tag'),u.Amp);\n";
print $mf "set(xedt,'String',str);\n";
print $mf "disp (str);\n";
print $mf "end\n\n";


print $mf "function BDFssp(h,e)\nu=get(h,'UserData');\nfprintf('%s\\n%s\\n',get(h,'Tag'),u.text);\nset(xedt,'String',get(h,'Tag'));\nend\n\n";
print $mf "\n\nend\n";
