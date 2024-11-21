function knots_new = HB_subdivide(knots,unacceptable,m)
%SUBDIVIDE 对于自适应B样条，判断哪些地方需要加密，并将相应区间的中点纳入节点序列
%   输入：knots：原有的节点序列
%   输入：unacceptable：出现容许以外误差需要进一步细分的区间编号
%   输出：knots_new：新的细分节点。在层次B样条中可以直接取原节点序列的所有中点

knots_new=knots;
insert_key=1;
for i=1:size(knots,1)-1
    if unacceptable(i,1)==1
        knots_new_before=knots_new(1:insert_key);
        knots_new_after=knots_new(insert_key+1:end);
        insertion=(knots(i)+knots(i+1))/2;
        knots_new=[knots_new_before;insertion;knots_new_after];
        insert_key=insert_key+2;
    else
        insert_key=insert_key+1;
    end
end

end

