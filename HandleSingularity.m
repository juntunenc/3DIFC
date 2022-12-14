%%
function new_Field = HandleSingularity(Field)

new_Field = Field;
new_Field(isnan(Field)==1) = 0;
new_Field(isinf(abs(Field))==1) = 0;