function num = bit2int_(bits, order)

num = sum(bits.*2.^[(order-1:-1:0)].',1);