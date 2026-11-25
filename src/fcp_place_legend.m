function info = fcp_place_legend(ax,lgd,options)
if nargin<3, options=struct(); end
defaults={'margin',[0.03,0.17];'clearance',0.06;'comfort',0.14;'step',0.004};
for k=1:size(defaults,1)
    if ~isfield(options,defaults{k,1}), options.(defaults{k,1})=defaults{k,2}; end
end
fig=ancestor(ax,'figure');
drawnow;
[axes_box,figure_size]=frame(ax,fig);
legend_size=normalized_size(lgd).*figure_size;
scale=[axes_box(3)*figure_size(1),axes_box(4)*figure_size(2)];
box_size=(legend_size+options.margin(2))./scale;
margin=options.margin(1)./scale;
if any(box_size+2*margin>=1)
    error('fcp:legend','The legend is larger than the axes.');
end
points=curve_points(ax,options.step);

names={'northeast','northwest','southeast','southwest','north','south','east','west'};
low=margin; high=1-margin-box_size; middle=(1-box_size)/2;
corners=[high(1),high(2);low(1),high(2);high(1),low(2);low(1),low(2); ...
    middle(1),high(2);middle(1),low(2);high(1),middle(2);low(1),middle(2)];
distance=zeros(size(corners,1),1);
for k=1:size(corners,1)
    distance(k)=box_distance(points,corners(k,:),box_size,scale);
end
info.raised=false; info.ylim=get(ax,'YLim');
choice=find(distance>=options.comfort,1);
if isempty(choice)
    [best,choice]=max(distance);
    if best<options.clearance, choice=[]; end
end
if isempty(choice)
    top=[1,2,5]; factor=Inf(size(top));
    room=1-box_size(2)-margin(2)-options.comfort/scale(2);
    if room<=0.2
        error('fcp:legend','The legend is too tall for the axes.');
    end
    reach=options.comfort/scale(1);
    for k=1:numel(top)
        left=corners(top(k),1)-reach; right=corners(top(k),1)+box_size(1)+reach;
        inside=points(:,1)>=left & points(:,1)<=right;
        factor(k)=max([points(inside,2);0])/room;
    end
    [factor,k]=min(factor); choice=top(k);
    factor=max(factor,1);
    limits=get(ax,'YLim');
    if strcmp(get(ax,'YScale'),'log')
        limits(2)=10^(log10(limits(1))+factor*(log10(limits(2))-log10(limits(1))));
    else
        limits(2)=limits(1)+factor*(limits(2)-limits(1));
    end
    set(ax,'YLim',limits);
    info.raised=true; info.ylim=limits;
    drawnow;
    points=curve_points(ax,options.step);
    distance(choice)=box_distance(points,corners(choice,:),box_size,scale);
end
set(lgd,'Location',names{choice});
drawnow;
info.location=names{choice};
info.distance=distance(choice);
info.box_in_axes=[corners(choice,:),box_size];
setappdata(lgd,'fcp_legend_placement',info);
end

function [axes_box,figure_size]=frame(ax,fig)
old=get(ax,'Units'); set(ax,'Units','normalized');
axes_box=get(ax,'Position'); set(ax,'Units',old);
old=get(fig,'Units'); set(fig,'Units','inches');
position=get(fig,'Position'); set(fig,'Units',old);
figure_size=position(3:4);
end

function value=normalized_size(lgd)
old=get(lgd,'Units'); set(lgd,'Units','normalized');
position=get(lgd,'Position'); set(lgd,'Units',old);
value=position(3:4);
end

function points=curve_points(ax,step)
lines=findobj(ax,'Type','line');
points=zeros(0,2);
for k=1:numel(lines)
    x=frame_coordinate(get(lines(k),'XData'),get(ax,'XLim'),get(ax,'XScale'));
    y=frame_coordinate(get(lines(k),'YData'),get(ax,'YLim'),get(ax,'YScale'));
    x=x(:); y=y(:);
    if numel(x)~=numel(y) || isempty(x), continue; end
    if strcmp(get(lines(k),'LineStyle'),'none') || numel(x)==1
        points=[points;x,y];
        continue
    end
    count=max(1,ceil(hypot(diff(x),diff(y))/step));
    count(~isfinite(count))=0;
    total=sum(count)+numel(x);
    px=zeros(total,1); py=zeros(total,1); at=0;
    for j=1:numel(x)-1
        if count(j)==0, continue; end
        t=(0:count(j)).'/count(j); ids=at+(1:numel(t));
        px(ids)=x(j)+t*(x(j+1)-x(j)); py(ids)=y(j)+t*(y(j+1)-y(j));
        at=at+numel(t);
    end
    points=[points;px(1:at),py(1:at)];
end
points=points(all(isfinite(points),2),:);
points=points(points(:,1)>=0 & points(:,1)<=1 & points(:,2)>=0,:);
end

function t=frame_coordinate(values,limits,axis_scale)
values=double(values);
if strcmp(axis_scale,'log')
    t=NaN(size(values)); positive=values>0;
    t(positive)=(log10(values(positive))-log10(limits(1)))/(log10(limits(2))-log10(limits(1)));
else
    t=(values-limits(1))/(limits(2)-limits(1));
end
end

function distance=box_distance(points,corner,box_size,scale)
points=points(points(:,2)<=1,:);
if isempty(points), distance=Inf; return; end
dx=max(max(corner(1)-points(:,1),points(:,1)-corner(1)-box_size(1)),0)*scale(1);
dy=max(max(corner(2)-points(:,2),points(:,2)-corner(2)-box_size(2)),0)*scale(2);
distance=min(hypot(dx,dy));
end
