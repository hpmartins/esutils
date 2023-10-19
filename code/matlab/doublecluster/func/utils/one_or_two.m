function which_one_or_two = one_or_two(one, two, q_one_or_two)
    if q_one_or_two == 1
        which_one_or_two = one;
    elseif q_one_or_two == 2
        which_one_or_two = two;
    else
        which_one_or_two = 0;
    end
end