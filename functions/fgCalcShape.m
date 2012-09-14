function [curvature torsion] = fgCalcShape(fg)
% Calculate curvature and torsion for each node on each fiber


for jj = 1:length(fg)
    numfibers=length(fg(jj).fibers);
    for ii = 1:numfibers
        [~,~,~,curvature{jj,ii},torsion{jj,ii}] = frenet2(fg(jj).fibers{ii}(1,:)',fg(jj).fibers{ii}(2,:)',fg(jj).fibers{ii}(3,:)');
    end
end
