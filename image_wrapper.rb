require 'json'

cmd = ARGV.size > 0 ? ARGV.join(" ") : "-a result_b1_md.bam -s -S -p -r chr1:8869816-8899900 -l -U -y 20"
bin = File.exist?("./target/debug/hgb") ? "./target/debug/hgb" : "./target/release/hgb"
READ_MAX=20

File.open("./dnd/reads.json", 'w') do |file|
  JSON.dump({"read_max"=>READ_MAX} , file)
end
# Ruler
`#{bin} -t12 vis -_ 10000 #{cmd} -o dnd/0.png`
# Each read
for i in 0..READ_MAX
    `#{bin} -t12 vis -_ #{i} -o dnd/#{i+1}.png #{cmd} -*`
end
# Coverage
`#{bin} -t12 vis #{cmd} -A -o dnd/#{READ_MAX+1}.png`
