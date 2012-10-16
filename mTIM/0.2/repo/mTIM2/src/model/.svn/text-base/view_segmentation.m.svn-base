function view_segmentation(label_seq, LABELS, y_anchor, color)

% view_segmentation(label_seq, LABELS, y_anchor, color)
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2009

% draw gray background
X = [length(label_seq)+0.5;
    length(label_seq)+0.5;
    0.5;
     0.5];
Y = [y_anchor+0.22;
     y_anchor-0.22;
     y_anchor-0.22;
     y_anchor+0.22];
fill(X, Y, [0.7 0.7 0.7]);

% draw blocks corresponding to non-intergenic and non-intronic probes
exons_max = find_blocks(label_seq~=LABELS.intron_W ...
                        & label_seq~=LABELS.intron_C ...
                        & label_seq~=LABELS.intergenic);
num_exons = size(exons_max, 2);
X = [exons_max(2,:)+0.5;
     exons_max(2,:)+0.5;
     exons_max(1,:)-0.5;
     exons_max(1,:)-0.5];
Y = [(y_anchor+0.2)*ones(1,num_exons);
     (y_anchor-0.2)*ones(1,num_exons);
     (y_anchor-0.2)*ones(1,num_exons);
     (y_anchor+0.2)*ones(1,num_exons)];
fill(X, Y, color + [0.3 0.5 0.5]);

% draw blocks corresponding to exons on Watson strand
exons_min_W = find_blocks(label_seq==LABELS.exon_W);
num_exons_W = size(exons_min_W, 2);
X = [exons_min_W(2,:)+0.5;
     exons_min_W(2,:)+0.5;
     exons_min_W(1,:)-0.5;
     exons_min_W(1,:)-0.5];
Y = [(y_anchor+0.2)*ones(1,num_exons_W);
     y_anchor*ones(1,num_exons_W);
     y_anchor*ones(1,num_exons_W);
     (y_anchor+0.2)*ones(1,num_exons_W)];
fill(X, Y, color);

% draw blocks corresponding to exons on Crick strand
exons_min_C = find_blocks(label_seq==LABELS.exon_C);
num_exons_C = size(exons_min_C, 2);
X = [exons_min_C(2,:)+0.5;
     exons_min_C(2,:)+0.5;
     exons_min_C(1,:)-0.5;
     exons_min_C(1,:)-0.5];
Y = [y_anchor*ones(1,num_exons_C);
     (y_anchor-0.2)*ones(1,num_exons_C);
     (y_anchor-0.2)*ones(1,num_exons_C);
     y_anchor*ones(1,num_exons_C)];
fill(X, Y, color);


% draw blocks corresponding to non-intergenic and non-exonic probes
introns_max = find_blocks(label_seq~=LABELS.exon_W ...
                          & label_seq~=LABELS.exon_C ...
                          & label_seq~=LABELS.intergenic);
num_introns = size(introns_max, 2);
X = [introns_max(2,:)+0.5;
     introns_max(2,:)+0.5;
     introns_max(1,:)-0.5;
     introns_max(1,:)-0.5];
Y = [(y_anchor+0.05)*ones(1,num_introns);
     (y_anchor-0.05)*ones(1,num_introns);
     (y_anchor-0.05)*ones(1,num_introns);
     (y_anchor+0.05)*ones(1,num_introns)];
fill(X, Y, color + [0.5 0.5 0.3]);

% draw blocks corresponding to introns on Watson strand
introns_min_W = find_blocks(label_seq==LABELS.intron_W);
num_introns_W = size(introns_min_W, 2);
X = [introns_min_W(2,:)+0.5;
     introns_min_W(2,:)+0.5;
     introns_min_W(1,:)-0.5;
     introns_min_W(1,:)-0.5];
Y = [(y_anchor+0.05)*ones(1,num_introns_W);
     y_anchor*ones(1,num_introns_W);
     y_anchor*ones(1,num_introns_W);
     (y_anchor+0.05)*ones(1,num_introns_W)];
fill(X, Y, [0.0 0.0 0.0]);

% draw blocks corresponding to introns on Crick strand
introns_min_C = find_blocks(label_seq==LABELS.intron_C);
num_introns_C = size(introns_min_C, 2);
X = [introns_min_C(2,:)+0.5;
     introns_min_C(2,:)+0.5;
     introns_min_C(1,:)-0.5;
     introns_min_C(1,:)-0.5];
Y = [y_anchor*ones(1,num_introns_C);
     (y_anchor-0.05)*ones(1,num_introns_C);
     (y_anchor-0.05)*ones(1,num_introns_C);
     y_anchor*ones(1,num_introns_C)];
fill(X, Y, [0.0 0.0 0.0]);
