function tet2_coords = general_tet_mapping(T1, T2, tet1_coords)

    %tet1 ===================

    a1x = T1(1,1); %vertice 1
    a1y = T1(1,2);
    a1z = T1(1,3);

    b1x = T1(2,1); %vertice 2
    b1y = T1(2,2);
    b1z = T1(2,3);

    c1x = T1(3,1); %vertice 3
    c1y = T1(3,2);
    c1z = T1(3,3);

    d1x = T1(4,1); %vertice 4
    d1y = T1(4,2);
    d1z = T1(4,3);

    %tet2 ====================

    a2x = T2(1,1); %vertice 1
    a2y = T2(1,2);
    a2z = T2(1,3);

    b2x = T2(2,1); %vertice 2
    b2y = T2(2,2);
    b2z = T2(2,3);

    c2x = T2(3,1); %vertice 3
    c2y = T2(3,2);
    c2z = T2(3,3);

    d2x = T2(4,1); %vertice 4
    d2y = T2(4,2);
    d2z = T2(4,3);

    D1 = [d1x; d1y; d1z];
    D2 = [d2x; d2y; d2z];

    A1 = [a1x b1x c1x; 
          a1y b1y c1y; 
          a1z b1z c1z] - D1;

    A2 = [a2x b2x c2x; 
          a2y b2y c2y; 
          a2z b2z c2z] - D2;
    
    tet2_coords = A2 * (A1 \ (tet1_coords - D1)) + D2;

end

