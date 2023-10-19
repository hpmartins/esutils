function barplot(x, y, width, color)
    line([x'; x'], [zeros(1,size(x,1)); y'], 'LineWidth', width, 'Color', color, 'LineSmoothing', 'on');
end