# Reed–Solomon error correction(BCH) (encoder, decoder)

# Intruduction 
Reed–Solomon codes are a group of error-correcting codes that were introduced by Irving S. Reed and Gustave Solomon in 1960. They have many applications, the most prominent of which include consumer technologies such as MiniDiscs, CDs, DVDs, Blu-ray discs, QR codes, data transmission technologies such as DSL and WiMAX, broadcast systems such as satellite communications, DVB and ATSC, storage systems such as RAID 6, coding schemes used by NASA missions.

# Used algorithms
- Binary Gauss-Jordan Elimination
- Gauss-Jordan Elimination Polynomials in Galois field
- Detrminant of Matrix n x n
- Transpose Matrix
- Zech's logarithm
- Adding polynomials 
- Multiplication polynomials
- Dividing polynomials
- Peterson–Gorenstein–Zierler decoder
- Chien search
- Cryptographic PRNGs (finding primitive polynomials)

# Instalation
```
yarn install
yarn start
```

# Encoder
- Finding primitive polynomial
- Finding elements of finite field
- Finding minimals polynomials
- Generator polynomial G(x)
- Encoding 
    - (m(x) * x^n-k) / G(x) then take Remainder R(x)
    - x^n-k + R(x)

# Decoder
- Syndorm -> (x^n-k + R(x)) / R(x)
- if Syndrom == 0 then end;
- if Syndrom != 0 then: 
    - make Matrix of Syndorms S
    - calculate determinat then det(S) 
    - loop if det(S) == 0 then decrease matrix rank abot 1 end;
    - if (det(S) != 0) then find result of equation 
- Result of equations are coefficients of Lambda Function
- Find Roots of Lambda functions 
- (Roots)^-1 are positions of wrong bit in code

