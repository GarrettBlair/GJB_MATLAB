function [neframe] =  quick_NE_reverse(frame_in, structuring_element)


    neframe = anisodiff2D(frame_in, 5, .25, .5, 1);
    bg = imclose(neframe, structuring_element);
    neframe = neframe-bg;

end