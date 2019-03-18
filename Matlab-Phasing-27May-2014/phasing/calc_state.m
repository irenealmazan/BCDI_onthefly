function flip=calc_state(trigger,qq)
%returns 1 or -1 if iteration number matches one of trigger elements
%if there is none returns 1 else returns -1
xx=fix(trigger)-fix(qq);
state=( xx == 0 );
state=sum(state);
if state == 1, flip = -1;, else flip = 1;,end;

end