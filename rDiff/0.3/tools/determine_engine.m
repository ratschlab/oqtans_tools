function [engine, environment, basedir, mccdir, license_type] = determine_engine()
% DETERMINES_ENGINE   Determines used interpreter.
%
%   [engine, environment, basedir, mccdir, license_type] = determine_engine() ;
% 
%   returns either 'matlab' or 'octave' in variable engine
%   and 'internal' or 'galaxy' in variable environment
%   basedir is the directory of the matlab or octave instance
%   mccdir is the directory of the matlab compiler (does not yet exist for octave)
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


global g_license_type ;

if isempty(g_license_type),
  g_license_type = 'academic' ;
end ;
license_type = g_license_type ;

lic=license ;
if ~isequal(lic, 'GNU General Public License'),
  engine='matlab' ;
else
  engine='octave' ;
end ;

environment='internal' ;
if isequal(whoami, 'galaxy'),
  environment='galaxy' ;
end ;

if isequal(engine, 'matlab') && isequal(environment,  'internal'),
  basedir = '' ;
elseif isequal(engine, 'matlab') && isequal(environment,  'galaxy'),
  basedir = '' ;
elseif isequal(engine, 'octave') && isequal(environment,  'internal'),
  basedir = '' ;
elseif isequal(engine, 'octave') && isequal(environment,  'galaxy'),
  basedir = '' ;
end ;

if isequal(environment,  'internal'),
  mccdir = '' ;
elseif isequal(environment,  'galaxy'),
  mccdir = '' ;
end ;

return
