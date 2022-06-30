function [neframe] =  quick_NE(frame_in, structuring_element)


    neframe = anisodiff2D(frame_in, 5, .25, .5, 1);
    bg = imopen(neframe, structuring_element);
    neframe = neframe-bg;

end