function X = GenLS(X,p)
ind = find(rand([size(X),1])<p);
X(ind) = rand(length(ind),1);
end