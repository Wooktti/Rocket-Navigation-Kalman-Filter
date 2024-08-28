function q_output = quaternion_product(q2,q1)

A = [q2(4), q2(3), -q2(2), q2(1);
    -q2(3), q2(4), q2(1),  q2(2);
    q2(2),  -q2(1), q2(4),  q2(3);
    -q2(1), -q2(2), -q2(3), q2(4)];
q_output = A * q1;
end