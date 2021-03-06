function [B, U, para,results] = prepare_onesite(A, direction,para,sitej,results)

[D1, D2, d] = size(A);
switch direction
    case 'lr'
        if para.parity=='n'
            A = permute(A, [3, 1, 2]);
            A = reshape(A, [d * D1, D2]);
            [B, S, U] = svd2(A);
            vNE = vonNeumannEntropy(S);
            sv = diag(S);
            DB = size(S, 1);
            para.D(sitej)=DB;
            B = reshape(B, [d, D1, DB]);
            B = permute(B, [2, 3, 1]);
            U = S * U;
        else
            A=permute(A,[3,1,2]);
            s=A(para.Anzilr{sitej});
            %s=nonzeros(A);
            sdim=length(s);
            assert(sdim==D1*D2*d/2);
            A_odd=reshape(s(1:sdim/2),[D1*d/2,D2/2]);
            A_even=reshape(s(sdim/2+1:end),[D1*d/2,D2/2]);
            [B_odd,S_odd,U_odd]=svd2(A_odd);
            [B_even,S_even,U_even]=svd2(A_even);
            
            vNE_odd=vonNeumannEntropy(S_odd);
            vNE_even=vonNeumannEntropy(S_even);
            vNE=vNE_odd+1i*vNE_even;
            DB_odd=size(S_odd,1);
            DB_even=size(S_even,1);
            DB=DB_odd+DB_even;
            sv_odd=diag(S_odd);
            sv_even=diag(S_even);
            sv=sv_odd+1i*sv_even;
            
            B=zeros(d,D1,DB);
            if D1~=1;
                for m=1:DB/2
                    submatsize=D1*d/4;
                    B(d/2+1:end,1:D1/2,m)=reshape(B_odd(1:submatsize,m),d/2,D1/2);
                    B(1:d/2,D1/2+1:end,m)=reshape(B_odd(submatsize+1:end,m),d/2,D1/2);
                    
                    B(1:d/2,1:D1/2,m+DB/2)=reshape(B_even(1:submatsize,m),d/2,D1/2);
                    B(d/2+1:end,D1/2+1:end,m+DB/2)=reshape(B_even(submatsize+1:end,m),d/2,D1/2);
                end
            else %Suppose the left most parity is even
                assert(DB==2);
                for m=1:DB/2
                    if para.parity=='e'
                        B(1:d/2,1,m)=B_odd(:,m);
                        B(d/2+1:end,1,m+DB/2)=B_even(:,m);
                    else
                        B(2,1,1)=B_odd;
                        B(1,1,2)=B_even;
                    end
                end
            end
            B=permute(B,[2,3,1]);
            U=blkdiag(S_odd*U_odd,S_even*U_even);

        end
    case 'rl'
        if para.parity=='n'
            A = permute(A, [1, 3, 2]);
            A = reshape(A, [D1, d * D2]);
            [U, S, B] = svd2(A);
            vNE = vonNeumannEntropy(S);
            sv = diag(S);
            DB = size(S, 1);
            B = reshape(B, [DB, d, D2]); B = permute(B, [1, 3, 2]);
            U = U * S;
        else
            A=permute(A,[2,3,1]);
            s=A(para.Anzirl{sitej});
            %s=nonzeros(A);
            sdim=length(s);
            assert(sdim==D1*D2*d/2);
            A_odd=reshape(s(1:sdim/2),[D2*d/2,D1/2]);
            A_even=reshape(s(sdim/2+1:end),[D2*d/2,D1/2]);
            [B_odd,S_odd,U_odd]=svd2(A_odd);
            [B_even,S_even,U_even]=svd2(A_even);
            
            vNE_odd=vonNeumannEntropy(S_odd);
            vNE_even=vonNeumannEntropy(S_even);
            vNE=vNE_odd+1i*vNE_even;
            DB_odd=size(S_odd,1);
            DB_even=size(S_even,1);
            DB=DB_odd+DB_even;
            sv_odd=diag(S_odd);
            sv_even=diag(S_even);
            sv=sv_odd+1i*sv_even;
            
            B=zeros(D2,d,DB);
            if D2~=1;
                for m=1:DB/2
                    submatsize=D2*d/4;
                    B(D2/2+1:end,1:d/2,m)=reshape(B_odd(1:submatsize,m),D2/2,d/2);
                    B(1:D2/2,d/2+1:end,m)=reshape(B_odd(submatsize+1:end,m),D2/2,d/2);
                    
                    B(1:D2/2,1:d/2,m+DB/2)=reshape(B_even(1:submatsize,m),D2/2,d/2);
                    B(D2/2+1:end,d/2+1:end,m+DB/2)=reshape(B_even(submatsize+1:end,m),D2/2,d/2);
                end
            else %Suppose the right most parity is even
                for m=1:DB/2
                    B(1,1:d/2,m)=B_odd(:,m);
                    B(1,d/2+1:end,m+DB/2)=B_even(:,m);
                end
            end
            B=permute(B,[3,1,2]);
            U=blkdiag(S_odd*U_odd,S_even*U_even);
            U=transpose(U);
            

        end
end
if nargin==5
    results.Amat_vNE(sitej) =vNE;
    if length(sv)==para.D(sitej) || (length(sv)==para.D(sitej)/2 && (para.parity~='n'))
        results.Amat_sv{sitej}=sv;
    end
end

end
