function overlap = waveformOverlap( wf1, wf2, a )
%WAVEFORMOVERLAP Computes the wavefrom overlap between two waveforms.
%
% DESCRIPTION:
%   Computes the wavefrom overlap between two waveforms. The overlap metric
%   was created for finding the similarity in two DIRSIG-generated waveform
%   signals. DIRSIG produces a simulated mean number of photons per time
%   bin. Waveform overlap is the intersection divided by the union of the
%   (1-a)*100% central regions (as defined by the Poisson distrubition)
%   from both waveforms.
%
% SYNTAX:
%   overlap = waveformOverlap(wf1, wf2)
%   overlap = waveformOverlap(wf1, wf2, alpha)
%
% INPUTS:
%   wf1: the first waveform
%   wf2: the second waveform
%   a: (optional) the alpha value
%
% OUTPUTS:
%   overlap: the overlap between the two waveforms.
%
% LICENSE:
%   The MIT License (MIT)
%
%   Copyright (c) 2012-2015 Rochester Institute of Technology
%
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to
%   permit persons to whom the Software is furnished to do so, subject to
%   the following conditions:
%
%   The above copyright notice and this permission notice shall be included
%   in all copies or substantial portions of the Software.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% AUTHOR(S):
%   Paul Romanczyk (par4249 at rit dot edu)
%
% PUBLIC REPOSITORY:
%   http://github.com/pavdpr/dissertationcode
%
% REFERENCES:
%   [1] Romanczyk, P., van Aardt, J., Cawse-Nicholson, K., Kelbe, D., 
%   McGlinchy, J., and Krause, K. (2013). Assessing the Impact of Broadleaf
%   Tree Structure on Airborne Full-waveform Small-footprint Lidar Signals
%   Through Simulation. Canadian Journal of Remote Sensing. 39 (S1), pp. 
%   S60?S72. DOI: 10.5589/m13-015.
%

    % check the inputs
    if ( nargin < 3 )
        a = 0.05;
    end
    s1 = size( wf1 );
    s2 = size( wf2 );
    if ( numel( s1 ) ~= numel( s2 ) )
        error( 'wf1 and wf2 must be the same size' );
    end
    if ( sum( abs( s1 - s2 ) ) ~= 0 )
        error( 'wf1 and wf2 must be the same size' );
    end
    
    % compute the limits of the waveform
    info1 = waveformLimits( wf1, a );
    info2 = waveformLimits( wf2, a );

    % find the upper and lower regions
    u = zeros( [ size( wf1 ), 2 ] );
    l = u;
    u( :, :, 1 ) = info1.u;
    u( :, :, 2 ) = info2.u;
    l( :, :, 1 ) = info1.l;
    l( :, :, 2 ) = info2.l;
    
    u = min( u, [], 3 );
    l = max( l, [], 3 );
    
    % find the overlap
    diff = u - l;
    
    % only count it if it is possitive
    diff( diff <= 0 ) = 0;
    
    % compute the areas
    aintersection = sum( diff, 2 );
    a1 = sum( info1.u - info1.l, 2 );
    a2 = sum( info2.u - info2.l, 2 );
    aunion = a1 + a2 - aintersection;
    
    % compute the overlap
    overlap = aintersection ./ aunion;
end

function info = waveformLimits( wf, a )
    % computes the upper and lower bounds on a waveform
    info.u = poissinv( 1 - ( a / 2 ), wf );
    info.l = poissinv( a / 2, wf );
end