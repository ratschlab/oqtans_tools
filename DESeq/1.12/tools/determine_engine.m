function engine = determine_engine()

lserve=license;
if ~isequal(lserve, 'GNU General Public License'),
    engine='matlab';
else
    engine='octave';
end;

return 
