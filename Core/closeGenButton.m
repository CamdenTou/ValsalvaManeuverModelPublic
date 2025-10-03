function [] = closeGenButton(a,b)
% function closeGenButton
% closes a current figure window
% Input: a handle to figure, b not used
  f = a.Parent;
  close(f)
end 
