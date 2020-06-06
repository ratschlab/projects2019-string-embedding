function X = gen_strings(siglen,N,L,mutate_rate,block_mute,max_blocks)

% Generate data
X = cell(N,1); % inputs
x_prev = randi(siglen,1,L); % random  T x D
for i=1:N
    % point mutation
    x = x_prev; 
    mutate_inds = rand(1,L)<mutate_rate;
    num_mut = sum(mutate_inds);
    x(mutate_inds) = randi(siglen,1,num_mut);
    if rand<block_mute
        l = length(x);
        num_blocks = randi(max_blocks)+1;
        l2 = ceil(l/num_blocks)*num_blocks;
        x(l+1:l2) = randi(siglen,1,l2-l);
     
        rb_perm = randperm(num_blocks);
        blen = round(length(x)/num_blocks);
        x2 = zeros(size(x));
        for bi = 1:num_blocks
            br = rb_perm(bi);
            x2((1:blen) + (bi-1)*blen) = x((1:blen) + (br-1)*blen);
        end
        x = x2;   
    end
    
    
    x_prev = x; % save previous 
    X{i} = x; % save sample i 
end
