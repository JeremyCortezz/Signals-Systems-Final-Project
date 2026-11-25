function h = shade(ax, iv, color, alpha)
   if isempty(iv), h = gobjects(0); return; end
   K = size(iv, 1);
   X = [iv(:, 1)'; iv(:, 2)'; iv(:, 2)'; iv(:, 1)'];
   Y = repmat([-1e9; -1e9; 1e9; 1e9], 1, K);
   h = patch(ax, X, Y, color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
   uistack(h, 'bottom');
end
