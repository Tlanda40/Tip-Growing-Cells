function fig2pretty(varargin)
%FIG2PRETTY   publication-quality MATLAB figure transformation & eps export
%   FIG2PRETTY will prettify the current figure handle (gcf), with
%   the option of exporting it to EPS format. Prior to calling this 
%   function, the user should prepare the figure axes titles, labels, legends, etc. 
% 
%   FIG2PRETTY(fig_handle) will prettify the axes associated with
%   fig_handle, but will not export the figure to EPS.
%   
%   FIG2PRETTY(fname) will prettify the current figure handle (gcf) and export it
%   to 'fname.eps' or 'fname' if already concatenated with the eps extension
%  
%   FIG2PRETTY(fig_handle,fname) will combine the functionality of the
%   above two cases.
   
%   This function was written based on advice from 'Loren on the art of MATLAB'
%         http://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/
%   and the MATLAB newsletter:
%         http://www.mathworks.com/company/newsletters/digest/june00/export/index.html
%   Alexandre Colavin | Summer 2012 | acolavin@stanford.edu


% Argument screen
if nargin == 1
    if ischar(varargin{1})
        mygcf = gcf;
        fname = varargin{1};
    else
        mygcf = varargin{1};
        fname = '';
    end
    
elseif nargin == 2
    if ischar(varargin{1})
        mygcf = varargin{2};
        fname = varargin{1};
    else
        mygcf = varargin{1};
        fname = varargin{2};
    end
elseif nargin > 2
    error('fig2pretty only accepts up to two arguments. Type "help fig2pretty" for details.');
else
    mygcf = gcf;
    fname = '';
end


% extension check on filename
if ~strcmp(fname,'')
    if ~strcmp(fname(end-3:end),'.eps');
        fname = [fname,'.eps'];
    end
end

allAxes = findall(mygcf,'type','axes');
myleg = legend;
allAxes(ismember(allAxes,myleg)) = []; % remove legend axes from modification
for mygcanum = 1:numel(allAxes)
    
    % General improvement of plot
    set(allAxes(mygcanum), ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'LineWidth'   , 1         , ...
        'FontName'    , 'Helvetica' ,...
        'FontSize'    , 16);
    
    % Update titles and labels, if any
    curtitle = get(get(allAxes(mygcanum),'Title'),'String');
    curxlabel = get(get(allAxes(mygcanum),'Xlabel'),'String');
    curylabel = get(get(allAxes(mygcanum),'Ylabel'),'String');
    if ~strcmp(curtitle,'')
        set(get(allAxes(mygcanum),'Title'), ...
            'FontName','Helvetica', ...
            'FontSize',16);
    end
    if ~strcmp(curxlabel,'')
        set(get(allAxes(mygcanum),'Xlabel'), ...
            'FontName','Helvetica', ...
            'FontSize',16);
    end
    if ~strcmp(curylabel,'')
        set(get(allAxes(mygcanum),'Ylabel'), ...
            'FontName','Helvetica', ...
            'FontSize',16);
    end
    
    % Go through all line objects and make them a little thicker
    allplots = findall(allAxes(mygcanum),'type','line');
    for plotnum = 1:numel(allplots)
        set(allplots(plotnum),'LineWidth',1);
    end
    
end
set(gcf,'color','white')

if ~strcmp(fname,'')
    print(mygcf,'-depsc2',fname);
    
end

legend off
end



    

