function knots_new = HB_subdivide(knots,unacceptable,m)
%SUBDIVIDE ��������ӦB�������ж���Щ�ط���Ҫ���ܣ�������Ӧ������е�����ڵ�����
%   ���룺knots��ԭ�еĽڵ�����
%   ���룺unacceptable�������������������Ҫ��һ��ϸ�ֵ�������
%   �����knots_new���µ�ϸ�ֽڵ㡣�ڲ��B�����п���ֱ��ȡԭ�ڵ����е������е�

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

