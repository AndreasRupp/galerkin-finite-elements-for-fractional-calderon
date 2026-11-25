function save_figure_pdf(fig,filename,varargin)
if isempty(fig), fig=gcf; end
if ~isscalar(fig) || ~isgraphics(fig,'figure')
    error('fcp:figure','A valid figure handle is required.');
end
p=inputParser;
p.addParameter('ContentType','vector',@(s)ischar(s) || isstring(s));
p.addParameter('Resolution',400,@(x)isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
p.addParameter('Append',false,@(x)islogical(x) && isscalar(x));
p.parse(varargin{:});
opts=p.Results;
filename=char(filename);
if numel(filename)<4 || ~strcmpi(filename(end-3:end),'.pdf')
    filename=[filename,'.pdf'];
end
folder=fileparts(filename);
if ~isempty(folder) && ~exist(folder,'dir'), mkdir(folder); end
if isempty(which('exportgraphics'))
    error('fcp:export','PDF export requires MATLAB R2020a or later.');
end
ax=findall(fig,'Type','axes');
ax=ax(strcmp(get(ax,{'Visible'}),'on'));
target=fig;
if numel(ax)==1, target=ax; end
args={'ContentType',char(opts.ContentType),'BackgroundColor','white','Append',opts.Append};
if ~strcmpi(opts.ContentType,'vector')
    args=[args,{'Resolution',opts.Resolution}];
end
if ~verLessThan('matlab','25.1')
    args=[args,{'Padding',6,'Units','points'}]; % A small margin surrounds the labels and legend.
end
drawnow;
exportgraphics(target,filename,args{:});
end
