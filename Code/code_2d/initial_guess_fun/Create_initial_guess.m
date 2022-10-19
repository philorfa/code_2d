function value_guess=Create_initial_guess(guess_number)
initial_bottom=[4.68 3.21 0.0881]';
initial_top=[17.3 19.4 0.535]';
k=guess_number-1;
value_guess=(initial_top-initial_bottom)/k*(0:k)+repmat(initial_bottom,1,guess_number);
end